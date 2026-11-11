import { cigarOps } from "../src/cigar.js";

export const DEFAULTS = {
  minVaf: 0.02,
  mapq: 20,
  baseq: 20,
  duplicates: false,
  sample: "",
};
export function prepareRead(read, options = DEFAULTS) {
  const o = { ...DEFAULTS, ...options },
    flags = read.flags || 0;
  if (
    flags & (4 | 256 | 512 | 2048) ||
    (!o.duplicates && flags & 1024) ||
    read.mq === 255 ||
    !(read.mq >= o.mapq)
  )
    return null;
  if (o.sample && read.sample !== o.sample) return null;
  if (!read.seq || read.seq === "*" || !read.cigar || read.cigar === "*")
    return null;
  let ops;
  try {
    ops = cigarOps(read.cigar);
  } catch {
    return null;
  }
  let p = read.start,
    q = 0;
  const runs = [];
  for (const [n, op] of ops) {
    runs.push({ p, q, n, op });
    if ("M=XDN".includes(op)) p += n;
    if ("M=XIS".includes(op)) q += n;
  }
  if (q !== read.seq.length || !ops.length || !Number.isSafeInteger(p))
    return null;
  return { ...read, seq: read.seq.toUpperCase(), runs, end: p, baseq: o.baseq };
}
function good(read, q) {
  const quality = read.qual?.[q];
  return (
    /^[ACGT]$/.test(read.seq[q] || "") &&
    Number.isFinite(quality) &&
    quality !== 255 &&
    quality >= read.baseq
  );
}
function baseAt(read, position) {
  for (const r of read.runs)
    if ("M=X".includes(r.op) && position >= r.p && position < r.p + r.n) {
      const q = r.q + position - r.p;
      return good(read, q) ? read.seq[q] : null;
    }
  return null;
}
// Reconstruct only the requested reference span. Skips, uncovered gaps and
// low-quality bases remain unknown, never reference. Both indel flanks must be
// observed by the SAME alignment; mates may bridge different variant sites.
function sliceAllele(read, start, end, boundaries) {
  let covered = start,
    allele = "";
  for (const r of read.runs) {
    if (r.op === "I") {
      if (boundaries ? r.p >= start && r.p <= end : r.p > start && r.p < end) {
        for (let q = r.q; q < r.q + r.n; q++) if (!good(read, q)) return null;
        allele += read.seq.slice(r.q, r.q + r.n);
      }
      continue;
    }
    if (!"M=XDN".includes(r.op)) continue;
    const s = Math.max(start, r.p),
      e = Math.min(end, r.p + r.n);
    if (e <= s) continue;
    if (s !== covered || r.op === "N") return null;
    if (r.op !== "D") {
      const q = r.q + s - r.p;
      for (let j = 0; j < e - s; j++) if (!good(read, q + j)) return null;
      allele += read.seq.slice(q, q + e - s);
    } else if (r.p < start || r.p + r.n > end) return null;
    covered = e;
  }
  return covered === end ? allele : null;
}
export function callVariant(read, variant) {
  if (read.chr !== variant.chr) return "?";
  const { pos, refLength: length, ref, alt } = variant,
    end = pos + length;
  const indel = ref === null || alt === null || length !== alt.length;
  if (indel && (!baseAt(read, pos - 1) || !baseAt(read, end))) return "?";
  const allele = sliceAllele(read, pos, end, indel);
  if (allele === null) return "?";
  if (alt === null) {
    if (allele === "") return "R";
    if (allele.length !== variant.copyLength) return "O";
    const source = sliceAllele(
      read,
      variant.copyStart,
      variant.copyStart + variant.copyLength,
      false,
    );
    return source === null ? "?" : source === allele ? "A" : "O";
  }
  if (allele === alt) return "A";
  if (ref !== null ? allele === ref : allele.length === length) return "R";
  return "O";
}
export function windowsForVariants(variants, contigs) {
  const lengths = new Map(contigs.map((c) => [c.refName, c.length]));
  const windows = variants
    .map((v) => ({
      chr: v.chr,
      start: Math.max(0, Math.min(v.pos - 1, v.copyStart ?? v.pos)),
      end: Math.min(
        lengths.get(v.chr) ?? Infinity,
        Math.max(
          v.pos + v.refLength + 1,
          (v.copyStart ?? 0) + (v.copyLength ?? 0),
        ),
      ),
    }))
    .sort((a, b) => a.chr.localeCompare(b.chr) || a.start - b.start);
  const merged = [];
  for (const w of windows) {
    const previous = merged.at(-1);
    if (
      previous &&
      previous.chr === w.chr &&
      w.start <= previous.end + 250 &&
      w.end - previous.start <= 200000
    )
      previous.end = Math.max(previous.end, w.end);
    else merged.push({ ...w });
  }
  return merged;
}
const lowerBound = (arr, value) => {
  let lo = 0,
    hi = arr.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid].pos < value) lo = mid + 1;
    else hi = mid;
  }
  return lo;
};
export function fragmentKey(r) {
  return JSON.stringify([
    r.sample || "",
    r.rg || "",
    r.name,
    r.flags & 1 ? "pair" : `single:${r.fileOffset ?? r.id}`,
  ]);
}
export function makeCollector(variants, options) {
  const blocked = new Set(),
    fragments = new Map(),
    byChr = new Map(),
    stats = { fetched: 0, filtered: 0, usable: 0 };
  let maxSpan = 1;
  for (const v of variants) {
    if (!byChr.has(v.chr)) byChr.set(v.chr, []);
    byChr.get(v.chr).push(v);
    maxSpan = Math.max(maxSpan, v.refLength);
  }
  for (const list of byChr.values()) list.sort((a, b) => a.pos - b.pos);
  return {
    stats,
    fragments,
    add(raw) {
      stats.fetched++;
      const rawKey = fragmentKey(raw);
      if (!options?.duplicates && raw.flags & 1024) {
        blocked.add(rawKey);
        const prior = fragments.get(rawKey);
        if (prior) {
          stats.usable -= prior.ends.length;
          stats.filtered += prior.ends.length;
          fragments.delete(rawKey);
        }
        stats.filtered++;
        return;
      }
      if (blocked.has(rawKey)) {
        stats.filtered++;
        return;
      }
      const read = prepareRead(raw, options);
      if (!read) {
        stats.filtered++;
        return;
      }
      stats.usable++;
      // QNAME is scoped by sample AND read group. Only primary paired segments
      // share a fragment; cross-contig ends cannot establish local cis phase.
      const key = fragmentKey(read);
      let f = fragments.get(key);
      if (!f) {
        f = {
          key,
          name: read.name,
          sample: read.sample || "",
          rg: read.rg || "",
          calls: new Map(),
          ends: [],
        };
        fragments.set(key, f);
      }
      if (f.ends.length >= 2 && read.flags & 1)
        throw Error(
          `More than two primary alignments for template ${read.name}; QNAME/read-group identity is ambiguous`,
        );
      if (
        read.flags & 1 &&
        f.ends.some((e) => e.endFlag === (read.flags & 192))
      )
        throw Error(
          `Repeated primary read end for ${read.name}; QNAME/read-group identity is ambiguous`,
        );
      if (read.flags & 1 && f.ends.length && read.mate >= 0) {
        const previous = f.ends[0];
        if (
          previous.start !== read.mate ||
          previous.chr !== read.mateChr ||
          (previous.mate >= 0 &&
            (previous.mate !== read.start || previous.mateChr !== read.chr))
        )
          throw Error(
            `Inconsistent paired-end mate coordinates for ${read.name}`,
          );
      }
      f.ends.push({
        mate: read.mate,
        mateChr: read.mateChr,
        chr: read.chr,
        start: read.start,
        end: read.end,
        endFlag: read.flags & 192,
        strand: read.flags & 16 ? "-" : "+",
        mq: read.mq,
        cigar: read.cigar,
      });
      const list = byChr.get(read.chr) || [];
      for (
        let i = lowerBound(list, read.start - maxSpan);
        i < list.length && list[i].pos <= read.end;
        i++
      ) {
        const v = list[i],
          call = callVariant(read, v);
        if (call === "?") continue;
        const old = f.calls.get(v.id);
        f.calls.set(v.id, old === undefined || old === call ? call : "X");
      }
    },
  };
}
export function summarize(
  variants,
  fragments,
  stats = {},
  pairBudget = 2000000,
) {
  const sites = variants.map((v) => ({
      id: v.id,
      alt: 0,
      ref: 0,
      other: 0,
      conflict: 0,
    })),
    pairs = new Map();
  let comparisons = 0;
  for (const f of fragments.values()) {
    const calls = [...f.calls].sort((a, b) => a[0] - b[0]);
    for (const [id, call] of calls)
      sites[id][{ A: "alt", R: "ref", O: "other", X: "conflict" }[call]]++;
    for (let i = 0; i < calls.length; i++)
      for (let j = i + 1; j < calls.length; j++) {
        const [a, x] = calls[i],
          [b, y] = calls[j];
        if (variants[a].chr !== variants[b].chr) continue;
        if (++comparisons > 20000000)
          throw Error(
            "More than 20 million fragment/site-pair comparisons. Analyze smaller batches; no partial result was reported.",
          );
        const key = `${a}:${b}`;
        let p = pairs.get(key);
        if (!p) {
          if (pairs.size >= pairBudget)
            throw Error(
              "More than two million observed variant pairs. Analyze smaller batches; no variants were silently omitted.",
            );
          p = {
            a,
            b,
            both: 0,
            leftOnly: 0,
            rightOnly: 0,
            neither: 0,
            other: 0,
            conflict: 0,
          };
          pairs.set(key, p);
        }
        if (x === "X" || y === "X") p.conflict++;
        else if (x === "O" || y === "O") p.other++;
        else if (x === "A" && y === "A") p.both++;
        else if (x === "A") p.leftOnly++;
        else if (y === "A") p.rightOnly++;
        else p.neither++;
      }
  }
  return {
    sites,
    pairs: [...pairs.values()],
    fragments: fragments.size,
    ...stats,
  };
}
export function pairMetrics(p) {
  const valid = p.both + p.leftOnly + p.rightOnly + p.neither;
  const atLeastOne = p.both + p.leftOnly + p.rightOnly;
  return {
    valid,
    atLeastOne,
    together: atLeastOne ? p.both / atLeastOne : null,
    separate: atLeastOne ? (p.leftOnly + p.rightOnly) / atLeastOne : null,
    bGivenA: p.both + p.leftOnly ? p.both / (p.both + p.leftOnly) : null,
    aGivenB: p.both + p.rightOnly ? p.both / (p.both + p.rightOnly) : null,
  };
}
