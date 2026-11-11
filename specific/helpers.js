import {
  prepareRead,
  makeCollector,
  summarize,
  fragmentKey,
} from "./phasing.js";

export function templateKey(r) {
  return fragmentKey(r);
}
export function groupMolecules(reads, options) {
  const templates = new Map(),
    molecules = new Map(),
    keys = new Map();
  const blocked = new Set(reads.filter((r) => r.flags & 1024).map(templateKey));
  for (const raw of reads) {
    if (blocked.has(templateKey(raw))) continue;
    const r = prepareRead(raw, { ...options, duplicates: false });
    if (!r) continue;
    const key = templateKey(r);
    if (!templates.has(key)) templates.set(key, []);
    templates.get(key).push(r);
  }
  for (const [template, ends] of templates) {
    const r = ends[0],
      umis = [...new Set(ends.map((e) => e.umi || "").filter(Boolean))];
    if (umis.length > 1)
      throw Error(`Conflicting molecule tags on paired template ${r.name}`);
    const coordinates =
      r.flags & 1 && !(r.flags & 8) && r.mate >= 0 && r.mateChr
        ? [
            [r.chr, r.start, !!(r.flags & 16)],
            [r.mateChr, r.mate, !!(r.flags & 32)],
          ].sort((a, b) => JSON.stringify(a).localeCompare(JSON.stringify(b)))
        : [[r.chr, r.start, r.end, !!(r.flags & 16)]];
    const molecule = JSON.stringify([
      r.sample,
      r.rg,
      coordinates,
      umis[0] || "",
    ]);
    keys.set(template, molecule);
    if (!molecules.has(molecule)) molecules.set(molecule, []);
    molecules.get(molecule).push(...ends);
  }
  return { molecules, keys };
}
// Discover biallelic aligned-base markers, counting templates once per position.
// Without a supplied genome, the modal observed base is the baseline (not a
// claim that it is the genomic reference allele).
export function discoverHelpers(reads, requested, options, limits = {}) {
  const minSupport = options.helperSupport ?? 5,
    minFraction = options.helperFraction ?? 0.02;
  const positions = new Map(),
    maxPositions = limits.maxPositions ?? 1000000;
  const excluded = new Map();
  for (const v of requested) {
    if (!excluded.has(v.chr)) excluded.set(v.chr, []);
    excluded.get(v.chr).push([v.pos - 1, v.pos + Math.max(v.refLength, 1) + 1]);
  }
  const { molecules } = groupMolecules(reads, options);
  for (const ends of molecules.values()) {
    const votes = new Map();
    for (const r of ends)
      for (const run of r.runs) {
        if (!"M=X".includes(run.op)) continue;
        for (let j = 0; j < run.n; j++) {
          const p = run.p + j,
            q = run.q + j,
            base = r.seq[q],
            quality = r.qual?.[q];
          if (
            !/^[ACGT]$/.test(base) ||
            !Number.isFinite(quality) ||
            quality === 255 ||
            quality < options.baseq
          )
            continue;
          const key = `${r.chr}:${p}`;
          let site = positions.get(key);
          if (!site) {
            if (positions.size >= maxPositions)
              throw Error(
                "Helper search exceeds one million covered bases. Reduce the helper search gap or disable helper SNPs.",
              );
            site = { chr: r.chr, pos: p, counts: [0, 0, 0, 0] };
            positions.set(key, site);
          }
          const old = votes.get(key);
          votes.set(key, old === undefined || old === base ? base : "?");
        }
      }
    for (const [key, base] of votes)
      if (base !== "?") positions.get(key).counts["ACGT".indexOf(base)]++;
  }
  const helpers = [];
  for (const site of positions.values()) {
    if (excluded.get(site.chr)?.some(([s, e]) => site.pos >= s && site.pos < e))
      continue;
    const alleles = site.counts
      .map((n, i) => ["ACGT"[i], n])
      .filter(([, n]) => n)
      .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
    if (alleles.length !== 2) continue;
    const [[ref, refCount], [alt, altCount]] = alleles,
      depth = refCount + altCount;
    if (
      altCount < minSupport ||
      refCount < minSupport ||
      altCount / depth <= minFraction
    )
      continue;
    helpers.push({
      id: requested.length + helpers.length,
      chr: site.chr,
      pos: site.pos,
      ref,
      alt,
      refLength: 1,
      helper: true,
      label: `${site.chr}:${site.pos + 1}${ref}>${alt}`,
      sources: [],
      altCount,
      depth,
      fraction: altCount / depth,
    });
    if (helpers.length > (limits.maxHelpers ?? 2000))
      throw Error(
        "More than 2,000 qualifying helper SNPs. Reduce helper search scope or raise support/fraction thresholds; no helpers were silently downsampled.",
      );
  }
  return helpers;
}
export function helperWindows(variants, contigs, maxGap = 50000) {
  const lengths = new Map(contigs.map((c) => [c.refName, c.length]));
  const sorted = [...variants].sort(
      (a, b) => a.chr.localeCompare(b.chr) || a.pos - b.pos,
    ),
    out = [];
  for (const v of sorted) {
    const w = {
      chr: v.chr,
      start: Math.max(0, Math.min(v.pos - 250, v.copyStart ?? v.pos)),
      end: Math.min(
        lengths.get(v.chr),
        Math.max(
          v.pos + v.refLength + 250,
          (v.copyStart ?? 0) + (v.copyLength ?? 0),
        ),
      ),
    };
    const last = out.at(-1);
    if (last && last.chr === w.chr && w.start - last.end <= maxGap)
      last.end = Math.max(last.end, w.end);
    else out.push(w);
  }
  if (out.reduce((n, w) => n + w.end - w.start, 0) > 2000000)
    throw Error(
      "Helper search exceeds two million reference bases. Reduce the maximum search gap or disable helper SNPs.",
    );
  // Keep each individual indexed request bounded.
  return out.flatMap((w) => {
    const chunks = [];
    for (let s = w.start; s < w.end; s += 200000)
      chunks.push({ chr: w.chr, start: s, end: Math.min(w.end, s + 200000) });
    return chunks;
  });
}
export function helperGraph(variants, pairs, options) {
  const graph = new Map(),
    support = options.helperSupport ?? 5,
    agreement = options.linkAgreement ?? 0.95;
  for (const p of pairs) {
    if (!variants[p.a].helper && !variants[p.b].helper) continue;
    if (p.other || p.conflict) continue;
    const same = p.both + p.neither,
      opposite = p.leftOnly + p.rightOnly,
      total = same + opposite;
    const sign = same >= opposite ? 1 : -1;
    const supported =
      sign === 1
        ? p.both >= support && p.neither >= support
        : p.leftOnly >= support && p.rightOnly >= support;
    if (!supported || Math.max(same, opposite) / total < agreement) continue;
    const edge = {
      a: p.a,
      b: p.b,
      sign,
      count: Math.max(same, opposite),
      total,
      agreement: Math.max(same, opposite) / total,
    };
    for (const [a, b] of [
      [p.a, p.b],
      [p.b, p.a],
    ]) {
      if (!graph.has(a)) graph.set(a, []);
      graph.get(a).push({ ...edge, to: b });
    }
  }
  for (const edges of graph.values())
    edges.sort((a, b) => b.agreement - a.agreement || b.total - a.total);
  return graph;
}
// A path describes compatible marker relationships, never an invented joint
// frequency. Detect contradictory cycles across the entire reachable block.
export function findHelperPath(graph, variants, a, b) {
  const queue = [a],
    visited = new Map([[a, { sign: 1, parent: null, edge: null }]]);
  let inconsistent = false;
  for (let i = 0; i < queue.length; i++) {
    const current = queue[i],
      state = visited.get(current);
    for (const edge of graph.get(current) || []) {
      const sign = state.sign * edge.sign,
        prior = visited.get(edge.to);
      if (prior) {
        if (prior.sign !== sign) inconsistent = true;
        continue;
      }
      visited.set(edge.to, { sign, parent: current, edge });
      queue.push(edge.to);
    }
  }
  if (!visited.has(b)) return { status: "unresolved", nodes: [], edges: [] };
  const nodes = [],
    edges = [];
  for (let n = b; n !== null; n = visited.get(n).parent) {
    nodes.push(n);
    if (visited.get(n).edge) edges.push(visited.get(n).edge);
  }
  nodes.reverse();
  edges.reverse();
  const helperIds = nodes.filter((id) => variants[id].helper);
  if (!helperIds.length) return { status: "unresolved", nodes: [], edges: [] };
  return {
    status: inconsistent ? "conflicting" : "linked",
    sign: visited.get(b).sign,
    nodes,
    edges,
    helperIds,
  };
}
export function analyzeReads(reads, variants, options) {
  const helpers = options.helpers
    ? discoverHelpers(reads, variants, options)
    : [];
  let all = [...variants, ...helpers],
    collector = makeCollector(all, options);
  for (const r of reads) collector.add(r);
  const minimum = options.minVaf ?? 0.02;
  if (!Number.isFinite(minimum) || minimum < 0 || minimum > 1)
    throw Error("Minimum VAF must be between 0 and 100%");
  const initial = summarize(all, collector.fragments, collector.stats),
    suppressedVariants = [];
  all = all
    .filter((v) => {
      const site = initial.sites[v.id],
        total = site.alt + site.ref + site.other;
      if (total && site.alt / total < minimum) {
        suppressedVariants.push({
          ...v,
          vaf: site.alt / total,
          altCount: site.alt,
          callable: total,
        });
        return false;
      }
      return true;
    })
    .map((v, id) => ({ ...v, id }));
  if (suppressedVariants.length) {
    collector = makeCollector(all, options);
    for (const r of reads) collector.add(r);
  }
  let helperPairs = [];
  if (helpers.length) {
    const { keys } = groupMolecules(reads, options),
      molecules = new Map();
    const helperCollector = makeCollector(all, {
      ...options,
      duplicates: false,
    });
    for (const r of reads) helperCollector.add(r);
    for (const f of helperCollector.fragments.values()) {
      if (!keys.has(f.key)) continue;
      const key = keys.get(f.key);
      if (!molecules.has(key)) molecules.set(key, { calls: new Map() });
      const calls = molecules.get(key).calls;
      for (const [id, call] of f.calls) {
        const old = calls.get(id);
        calls.set(id, old === undefined || old === call ? call : "X");
      }
    }
    helperPairs = summarize(all, molecules).pairs;
  }
  return {
    ...summarize(all, collector.fragments, collector.stats),
    helperPairs,
    variants: all,
    requestedCount: all.filter((v) => !v.helper).length,
    helperCount: all.filter((v) => v.helper).length,
    suppressedVariants,
    collector,
  };
}
