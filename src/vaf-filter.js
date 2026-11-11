import {
  prepareRecord,
  alignmentRuns,
  hasSequence,
} from "./alignment-geometry.js";
const key = (r) => JSON.stringify([r.sample || "", r.rg || "", r.name]);
const add = (m, p, v) => {
  if (!m.has(p)) m.set(p, new Set());
  m.get(p).add(v);
};
function eventKey(p, ref, alt, reference, origin) {
  let allele = ref || alt;
  if (reference)
    while (p > origin && reference[p - origin - 1] === allele.at(-1)) {
      p--;
      allele = allele.at(-1) + allele.slice(0, -1);
    }
  return `${p}:${ref ? allele : ""}>${alt ? allele : ""}`;
}
// Suppression removes observations (unknown), never rewrites ALT into REF.
// Counts are unique physical fragments; contradictory observations do not vote.
export function applyVafFilter(reads, region, reference, opts = {}) {
  const minimum = opts.minVaf ?? 0.02;
  if (!Number.isFinite(minimum) || minimum < 0 || minimum > 1)
    throw Error("Minimum VAF must be between 0 and 100%");
  for (const r of reads) {
    r.vafHiddenBases = [];
    r.vafHiddenInsertions = [];
    r.vafHiddenDeletions = [];
    r.vafAssemblyMask = [];
  }
  if (!minimum) return { suppressed: 0, minimum, hiddenVariants: [] };
  const fragments = new Map(),
    events = new Map(),
    eligible = opts.eligible;
  for (const r of reads) {
    prepareRecord(r);
    if (
      !hasSequence(r) ||
      !r.alignmentUsable ||
      r.flags & (4 | 256 | 512) ||
      r.mq === 255 ||
      r.mq < (opts.mapq ?? 20) ||
      (opts.sample && r.sample !== opts.sample) ||
      (r.chr && region.chr && r.chr !== region.chr) ||
      (!(opts.includeDuplicates ?? opts.duplicates) && r.flags & 1024) ||
      (eligible && !eligible(r))
    )
      continue;
    const k = key(r);
    if (!fragments.has(k))
      fragments.set(k, {
        bases: new Map(),
        events: new Set(),
        records: [],
      });
    const combined = fragments.get(k),
      f = { bases: new Map(), ins: new Map(), events: new Set(), spans: [] },
      good = (q) =>
        q >= 0 &&
        q < r.seq.length &&
        (!r.qual?.length ||
          r.qual[q] === 255 ||
          r.qual[q] >= (opts.baseq ?? 20));
    for (const run of alignmentRuns(r)) {
      const { op, length: n, refStart: p, queryStart: q } = run;
      if ("M=X".includes(op)) {
        f.spans.push([p, p + n]);
        for (let j = 0; j < n; j++) {
          const pos = p + j;
          if (
            pos < region.start ||
            pos >= region.end ||
            !good(q + j) ||
            !/[ACGT]/.test(r.seq[q + j])
          )
            continue;
          add(f.bases, pos, r.seq[q + j]);
          if (j && good(q + j - 1)) add(f.ins, pos, "");
        }
      } else if (
        (op === "I" || op === "D") &&
        p >= region.start &&
        p < region.end
      ) {
        const valid =
          op === "I"
            ? Array.from({ length: n }, (_, j) => good(q + j)).every(Boolean)
            : good(q - 1) && good(q);
        if (!valid) continue;
        const ref =
            op === "D"
              ? reference?.slice(p - region.start, p - region.start + n) ||
                "-".repeat(n)
              : "",
          alt = op === "I" ? r.seq.slice(q, q + n) : "",
          id = eventKey(p, ref, alt, reference, region.start);
        if (!events.has(id))
          events.set(id, {
            id,
            pos: +id.split(":")[0],
            op,
            length: n,
            anchors: new Set(),
          });
        events.get(id).anchors.add(p);
        f.events.add(id);
        if (op === "I") add(f.ins, +id.split(":")[0], id.split(">")[1]);
        else
          for (let pos = p; pos < Math.min(p + n, region.end); pos++)
            add(f.bases, pos, "-");
      }
    }
    // Event presence takes precedence over an equivalent shifted empty
    // boundary within this alignment; disagreement between mates stays unknown.
    for (const id of f.events) {
      const e = events.get(id);
      if (e.op === "I") {
        const p = e.pos;
        if (f.ins.has(p)) f.ins.get(p).delete("");
      }
    }
    combined.records.push(f);
    for (const [p, bs] of f.bases)
      for (const b of bs) add(combined.bases, p, b);
    for (const id of f.events) combined.events.add(id);
  }
  const counts = new Map();
  for (const f of fragments.values())
    for (const [p, bs] of f.bases) {
      if (!counts.has(p)) counts.set(p, { counts: new Map(), total: 0 });
      const c = counts.get(p);
      for (const b of bs) if (!c.counts.has(b)) c.counts.set(b, 0);
      if (bs.size === 1) {
        const b = [...bs][0];
        c.total++;
        c.counts.set(b, c.counts.get(b) + 1);
      }
    }
  // Visible low-quality/duplicate observations can have zero eligible support.
  // Include their allele identities so they cannot survive as raw display-only
  // clashes when other fragments establish a callable denominator.
  for (const r of reads) {
    if (!hasSequence(r)) continue;
    for (const run of alignmentRuns(r)) {
      const { op, length: n, refStart: p, queryStart: q } = run;
      if ("M=X".includes(op)) {
        for (let j = 0; j < n; j++) {
          const c = counts.get(p + j),
            b = r.seq[q + j];
          if (c && /^[ACGT]$/.test(b) && !c.counts.has(b)) c.counts.set(b, 0);
        }
      } else if (
        (op === "I" || op === "D") &&
        p >= region.start &&
        p < region.end
      ) {
        const ref =
          op === "D"
            ? reference?.slice(p - region.start, p - region.start + n) ||
              "-".repeat(n)
            : "";
        const alt = op === "I" ? r.seq.slice(q, q + n) : "";
        const id = eventKey(p, ref, alt, reference, region.start);
        if (!events.has(id))
          events.set(id, {
            id,
            pos: +id.split(":")[0],
            op,
            length: n,
            anchors: new Set(),
          });
        events.get(id).anchors.add(p);
      }
    }
  }
  const hidden = new Map();
  for (const [p, c] of counts) {
    const ref =
      reference?.[p - region.start] ||
      [...c.counts]
        .filter(([b]) => b !== "-")
        .sort((a, b) => b[1] - a[1])[0]?.[0];
    for (const [b, n] of c.counts)
      if (b !== "-" && b !== ref && c.total && n / c.total < minimum) {
        if (!hidden.has(p)) hidden.set(p, new Set());
        hidden.get(p).add(b);
      }
  }
  const hiddenEvents = new Set();
  for (const e of events.values()) {
    let total = 0,
      alt = 0;
    for (const f of fragments.values()) {
      const votes = new Set();
      for (const r of f.records) {
        if (r.events.has(e.id)) {
          votes.add("ALT");
          continue;
        }
        // A repeat-normalized event can move beyond a read's edge. Any
        // observed equivalent CIGAR anchor is valid for testing absence.
        const callable = [e.pos, ...e.anchors].some((p) => {
          const end = p + e.length;
          return e.op === "I"
            ? r.ins.get(p)?.size === 1
            : r.spans.some(([s, t]) => s < p && t > end) &&
                r.bases.get(p - 1)?.size === 1 &&
                r.bases.get(end)?.size === 1;
        });
        if (callable) votes.add("REF");
      }
      // Contradictory mates are unknown, never two independent votes.
      if (votes.size !== 1) continue;
      total++;
      if (votes.has("ALT")) alt++;
    }
    if (total && alt / total < minimum) hiddenEvents.add(e.id);
  }
  for (const r of reads) {
    if (!hasSequence(r)) continue;
    for (const run of alignmentRuns(r)) {
      const { op, length: n, refStart: p, queryStart: q } = run;
      if ("M=X".includes(op)) {
        for (let j = 0; j < n; j++)
          if (hidden.get(p + j)?.has(r.seq[q + j]))
            r.vafHiddenBases.push(q + j);
      } else if (
        (op === "I" || op === "D") &&
        p >= region.start &&
        p < region.end
      ) {
        const ref =
            op === "D"
              ? reference?.slice(p - region.start, p - region.start + n) ||
                "-".repeat(n)
              : "",
          alt = op === "I" ? r.seq.slice(q, q + n) : "";
        if (hiddenEvents.has(eventKey(p, ref, alt, reference, region.start))) {
          (op === "I" ? r.vafHiddenInsertions : r.vafHiddenDeletions).push(p);
          for (
            let i = Math.max(0, q - 61);
            i < Math.min(r.seq.length, q + (op === "I" ? n : 0) + 61);
            i++
          )
            r.vafAssemblyMask.push(i);
        }
      }
    }
    r.vafAssemblyMask.push(...r.vafHiddenBases);
  }
  return {
    minimum,
    hiddenVariants: [...hidden]
      .flatMap(([p, bs]) =>
        [...bs].map((b) => `${p}:${reference?.[p - region.start] || ""}>${b}`),
      )
      .concat([...hiddenEvents]),
    suppressed:
      [...hidden.values()].reduce((n, s) => n + s.size, 0) + hiddenEvents.size,
  };
}
export const hiddenBase = (r, q) => r.vafHiddenBases?.includes(q) || false;
export const hiddenEvent = (r, op, pos) =>
  (op === "I" ? r.vafHiddenInsertions : r.vafHiddenDeletions)?.includes(pos) ||
  false;
