import { applyVafFilter, hiddenBase } from "../src/vaf-filter.js";
import { prepareRecord, alignmentRuns } from "../src/alignment-geometry.js";
import { structuralEvidence } from "../src/structural.js";

export const reverseComplement = (s) =>
  [...s]
    .reverse()
    .map((b) => ({ A: "T", C: "G", G: "C", T: "A" })[b] || "N")
    .join("");
const fragmentKey = (r) => JSON.stringify([r.sample || "", r.rg || "", r.name]);
const intersects = (a, b) => {
  for (const x of a) if (b.has(x)) return true;
  return false;
};
const compatible = (a, b) => {
  for (const [k, v] of a) if (b.has(k) && b.get(k) !== v) return false;
  return true;
};

// A strand-symmetric de Bruijn graph. Reverse-complement handles are kept
// throughout construction and folded only for display, never counted twice.
export function assemble(reads, region, reference = null, settings = {}) {
  const opts = {
    minVaf: 0.02,
    k: 31,
    support: 3,
    mapq: 20,
    baseq: 20,
    duplicates: false,
    ...settings,
  };
  if (
    !Number.isInteger(opts.k) ||
    opts.k < 15 ||
    opts.k > 61 ||
    opts.k % 2 !== 1
  )
    throw Error("Choose an odd k between 15 and 61");
  const vafFilter = applyVafFilter(reads, region, reference, opts);
  const { k } = opts,
    nodes = [],
    index = new Map(),
    edges = new Map(),
    fragments = [],
    fIndex = new Map(),
    threads = [],
    issues = {};
  const node = (key) => {
    if (!index.has(key)) {
      if (nodes.length >= 200000)
        throw Error(
          "Assembly exceeds 200,000 oriented k-mers. Narrow the interval or increase k. No reads were downsampled.",
        );
      index.set(key, nodes.length);
      nodes.push({
        id: nodes.length,
        key,
        out: new Set(),
        in: new Set(),
        fs: new Set(),
        positions: [],
        refs: [],
      });
    }
    return index.get(key);
  };
  const addEdge = (a, b, f) => {
    nodes[a].out.add(b);
    nodes[b].in.add(a);
    const key = a + ":" + b;
    if (!edges.has(key)) {
      if (edges.size >= 400000)
        throw Error(
          "Assembly exceeds 400,000 oriented edges. Narrow the interval.",
        );
      edges.set(key, new Set());
    }
    edges.get(key).add(f);
  };
  let bases = 0;
  for (const r of reads) {
    prepareRecord(r);
    const outsideTarget =
      (r.chr && r.chr !== region.chr) ||
      r.start >= region.end ||
      r.end <= region.start;
    const reason =
      r.displayIssue ||
      (outsideTarget
        ? "outside-target"
        : r.flags & 256
          ? "secondary-alignment"
          : r.flags & 512
            ? "QC-failed"
            : !opts.duplicates && r.flags & 1024
              ? "duplicate-flag"
              : r.mq === 255 || r.mq < opts.mapq
                ? "mapping-quality"
                : null);
    r.exclusion = reason;
    if (reason) {
      issues[reason] = (issues[reason] || 0) + 1;
      continue;
    }
    bases += r.seq.length;
    if (bases > 2000000)
      throw Error(
        "More than two million usable read bases. Narrow the interval; no reads were downsampled.",
      );
    const key = fragmentKey(r);
    if (!fIndex.has(key)) {
      fIndex.set(key, fragments.length);
      fragments.push({
        id: fragments.length,
        key,
        name: r.name,
        reads: [],
        nodes: new Set(),
        calls: new Map(),
        conflict: false,
        features: new Map(),
      });
    }
    const f = fragments[fIndex.get(key)];
    f.reads.push(r.id);
    const coords = new Int32Array(r.seq.length).fill(-1),
      inScope = new Uint8Array(r.seq.length);
    for (const run of alignmentRuns(r))
      if ("M=X".includes(run.op)) {
        for (let j = 0; j < run.length; j++) {
          const q = run.queryStart + j,
            p = run.refStart + j,
            b = r.seq[q];
          coords[q] = p;
          if (p >= region.start && p < region.end) inScope[q] = 1;
          if (
            hiddenBase(r, q) ||
            p < region.start ||
            p >= region.end ||
            !"ACGT".includes(b) ||
            (r.qual?.length && r.qual[q] < opts.baseq)
          )
            continue;
          // This guard includes singleton SNPs, not just graph bubbles or sites
          // that pass a support threshold. One contradictory call is a veto.
          if (f.calls.has(p) && f.calls.get(p) !== b) f.conflict = true;
          f.calls.set(p, b);
        }
      } else if (
        (run.op === "I" || run.op === "S") &&
        run.refStart >= region.start &&
        run.refStart < region.end
      ) {
        inScope.fill(1, run.queryStart, run.queryStart + run.length);
      }
    let valid = 0,
      prev = null,
      path = [];
    for (let q = 0; q < r.seq.length; q++) {
      const b = r.seq[q];
      valid =
        !r.vafAssemblyMask?.includes(q) &&
        inScope[q] &&
        "ACGT".includes(b) &&
        (!r.qual?.length || r.qual[q] >= opts.baseq)
          ? valid + 1
          : 0;
      if (valid < k) {
        prev = null;
        if (path.length) {
          threads.push({ read: r.id, fragment: f.id, path });
          path = [];
        }
        continue;
      }
      const key = r.seq.slice(q - k + 1, q + 1),
        a = node(key),
        rev = node(reverseComplement(key));
      nodes[a].fs.add(f.id);
      nodes[rev].fs.add(f.id);
      f.nodes.add(a);
      f.nodes.add(rev);
      const p = coords[q - k + 1];
      if (p >= 0)
        for (const n of [a, rev])
          if (nodes[n].positions.length < 64) nodes[n].positions.push(p);
      if (prev !== null) {
        addEdge(prev, a, f.id);
        addEdge(rev, index.get(reverseComplement(nodes[prev].key)), f.id);
      }
      path.push(a);
      prev = a;
    }
    if (path.length) threads.push({ read: r.id, fragment: f.id, path });
  }
  // Preserve both ends (including excluded/unavailable mates) in the evidence
  // inspector without letting off-target mate sequence seed remote assembly.
  for (const f of fragments) f.reads = [];
  for (const r of reads) {
    const id = fIndex.get(fragmentKey(r));
    if (id !== undefined) fragments[id].reads.push(r.id);
  }
  if (!nodes.length && !vafFilter.suppressed)
    throw Error(
      "No quality-eligible sequence spans k bases. Lower k or quality thresholds, or choose another interval. All input records are retained.",
    );
  if (reference)
    for (let p = 0; p <= reference.length - k; p++) {
      const key = reference.slice(p, p + k),
        rev = reverseComplement(key);
      if (index.has(key)) nodes[index.get(key)].refs.push(region.start + p);
      if (index.has(rev)) nodes[index.get(rev)].refs.push(region.start + p);
    }
  const rc = nodes.map((n) => index.get(reverseComplement(n.key)));
  // Compact maximal nonbranching chains. Entire cycles remain explicit loops.
  const visited = new Set(),
    unitigs = [],
    loc = new Map(),
    bySequence = new Map();
  function compact(start) {
    const path = [];
    let a = start;
    while (!visited.has(a)) {
      visited.add(a);
      path.push(a);
      if (nodes[a].out.size !== 1) break;
      const b = [...nodes[a].out][0];
      if (nodes[b].in.size !== 1 || visited.has(b)) break;
      a = b;
    }
    const seq =
        nodes[path[0]].key +
        path
          .slice(1)
          .map((i) => nodes[i].key.at(-1))
          .join(""),
      rev = reverseComplement(seq),
      forward = seq <= rev,
      canonical = forward ? seq : rev;
    let u = bySequence.get(canonical);
    if (!u) {
      const ps = path.flatMap((i) => nodes[i].positions).sort((a, b) => a - b),
        fs = new Set(path.flatMap((i) => [...nodes[i].fs]));
      u = {
        id: unitigs.length,
        seq: canonical,
        length: canonical.length,
        support: path.reduce(
          (min, i) => Math.min(min, nodes[i].fs.size),
          Infinity,
        ),
        fragments: [...fs],
        position: ps.length ? ps[Math.floor(ps.length / 2)] : region.start,
        reference: reference
          ? path.filter((i) => nodes[i].refs.length).length / path.length
          : null,
        cyclic: nodes[path.at(-1)].out.has(path[0]),
        members: path.length,
      };
      bySequence.set(canonical, u);
      unitigs.push(u);
    }
    path.forEach((id, i) =>
      loc.set(id, {
        unitig: u.id,
        orientation: forward ? "+" : "-",
        offset: forward ? i : path.length - 1 - i,
      }),
    );
  }
  for (const n of nodes)
    if (
      !visited.has(n.id) &&
      (n.in.size !== 1 || nodes[[...n.in][0]].out.size !== 1)
    )
      compact(n.id);
  for (const n of nodes) if (!visited.has(n.id)) compact(n.id);
  const links = [],
    linkIndex = new Set();
  for (const [key, fs] of edges) {
    const [a, b] = key.split(":").map(Number),
      x = loc.get(a),
      y = loc.get(b);
    if (x.unitig === y.unitig && !unitigs[x.unitig].cyclic) continue;
    const code = `${x.unitig}${x.orientation}>${y.unitig}${y.orientation}`,
      twin = `${y.unitig}${y.orientation === "+" ? "-" : "+"}>${x.unitig}${x.orientation === "+" ? "-" : "+"}`,
      canonical = code < twin ? code : twin;
    if (linkIndex.has(canonical)) continue;
    linkIndex.add(canonical);
    links.push({
      from: x.unitig,
      to: y.unitig,
      fromOrientation: x.orientation,
      toOrientation: y.orientation,
      support: fs.size,
    });
  }
  // Only observed read traversals define bubble alleles. A graph walk alone
  // never creates a haplotype. Searches are deliberately bounded; complex
  // repeats and nested branches stay visible as unresolved topology.
  const bubbles = [],
    bubbleKeys = new Set(),
    suppressedGraphNodes = new Set();
  let complex = 0;
  for (const fork of nodes) {
    if (fork.out.size < 2 || fork.out.size > 4 || fork.fs.size < opts.support)
      continue;
    if (bubbles.length >= 128) {
      complex++;
      continue;
    }
    const reaches = [...fork.out].map((first) => {
      const distances = new Map([[first, 1]]),
        queue = [first];
      for (let z = 0; z < queue.length && z < 512; z++) {
        const a = queue[z],
          d = distances.get(a);
        if (d >= Math.max(2 * k + 20, 100)) continue;
        for (const b of nodes[a].out)
          if (b !== fork.id && !distances.has(b)) {
            distances.set(b, d + 1);
            queue.push(b);
          }
      }
      return distances;
    });
    const candidates = [...reaches[0]]
      .filter(([a]) => reaches.every((m) => m.has(a)))
      .sort(
        (a, b) =>
          Math.max(...reaches.map((m) => m.get(a[0]))) -
          Math.max(...reaches.map((m) => m.get(b[0]))),
      );
    if (!candidates.length) {
      complex++;
      continue;
    }
    const join = candidates[0][0],
      code = fork.id + ":" + join,
      twin = rc[join] + ":" + rc[fork.id],
      canonical = code < twin ? code : twin;
    if (bubbleKeys.has(canonical)) continue;
    bubbleKeys.add(canonical);
    const alleles = [],
      alleleIndex = new Map();
    for (const t of threads) {
      for (const path of [t.path, t.path.map((i) => rc[i]).reverse()]) {
        const i = path.indexOf(fork.id);
        if (i < 0) continue;
        const j = path.indexOf(join, i + 1);
        if (j < 0 || j - i > 2 * k + 120) continue;
        const ids = path.slice(i, j + 1),
          seq =
            fork.key +
            ids
              .slice(1)
              .map((a) => nodes[a].key.at(-1))
              .join("");
        if (!alleleIndex.has(seq)) {
          alleleIndex.set(seq, alleles.length);
          alleles.push({ id: alleles.length, seq, nodes: ids, fs: new Set() });
        }
        alleles[alleleIndex.get(seq)].fs.add(t.fragment);
      }
    }
    if (alleles.length < 2) continue;
    for (const a of alleles)
      a.diagnostic = new Set(
        a.nodes.filter(
          (n) => !alleles.some((b) => a !== b && b.nodes.includes(n)),
        ),
      );
    const bubble = {
      id: bubbles.length,
      fork: fork.id,
      join,
      alleles,
      position: unitigs[loc.get(fork.id).unitig].position,
      anchored: false,
    };
    // Reference annotation needs unique, forward oriented anchors and exact
    // intervening sequence. Reference k-mers never add graph nodes or edges.
    if (reference) {
      const positions = (key) => {
        const ps = [];
        let p = -1;
        while ((p = reference.indexOf(key, p + 1)) >= 0) {
          ps.push(p);
          if (ps.length > 1) break;
        }
        return ps;
      };
      const left = positions(fork.key),
        right = positions(nodes[join].key);
      if (
        left.length === 1 &&
        right.length === 1 &&
        right[0] > left[0] &&
        right[0] - left[0] < 2000
      ) {
        bubble.anchored = true;
        bubble.position = region.start + left[0] + k;
        bubble.referenceSequence = reference.slice(left[0], right[0] + k);
      } else {
        const left = positions(reverseComplement(nodes[join].key)),
          right = positions(reverseComplement(fork.key));
        if (
          left.length === 1 &&
          right.length === 1 &&
          right[0] > left[0] &&
          right[0] - left[0] < 2000
        ) {
          bubble.anchored = true;
          bubble.position = region.start + left[0] + k;
          bubble.referenceSequence = reverseComplement(
            reference.slice(left[0], right[0] + k),
          );
        }
      }
    }
    for (const a of alleles)
      a.isReference = bubble.referenceSequence
        ? a.seq === bubble.referenceSequence
        : null;
    const callable = new Set(alleles.flatMap((a) => [...a.fs])).size;
    bubble.alleles = alleles.filter(
      (a) => a.isReference || !callable || a.fs.size / callable >= opts.minVaf,
    );
    for (const a of alleles)
      if (!bubble.alleles.includes(a))
        for (const n of a.diagnostic) suppressedGraphNodes.add(n);
    bubble.alleles.forEach((a, i) => {
      a.id = i;
    });
    if (bubble.alleles.length >= 2) bubbles.push(bubble);
  }
  // Sequence-varying alignment positions provide a strict contradiction veto
  // even for singleton/error branches too weak or complex to form a bubble.
  const variable = new Map();
  for (const f of fragments)
    for (const [p, b] of f.calls) {
      if (!variable.has(p)) variable.set(p, new Set());
      variable.get(p).add(b);
    }
  for (const [p, bs] of variable) if (bs.size < 2) variable.delete(p);
  for (const f of fragments) {
    f.guard = new Map([...f.calls].filter(([p]) => variable.has(p)));
    for (const b of bubbles) {
      const hits = b.alleles.filter(
        (a) => a.fs.has(f.id) || intersects(f.nodes, a.diagnostic),
      );
      if (hits.length === 1) f.features.set(b.id, hits[0].id);
      else if (hits.length > 1) f.conflict = true;
    }
  }
  // Link features only through a physical fragment. Disconnected loci never
  // imply whole-region clone combinations. Partial evidence can remain shared.
  const parent = bubbles.map((_, i) => i),
    find = (i) => (parent[i] === i ? i : (parent[i] = find(parent[i])));
  // Inspired by uncertainty cuts in WhatsHap polyphase: one bridging read is
  // not sufficient to join otherwise disconnected phase blocks. This is a
  // simple fragment-count cut, not their fixed-ploidy threading optimizer.
  const phaseLinkMap = new Map();
  for (const f of fragments) {
    if (f.conflict) continue;
    const ids = [...f.features.keys()].sort((a, b) => a - b);
    for (let a = 0; a < ids.length; a++)
      for (let b = a + 1; b < ids.length; b++) {
        const key = ids[a] + ":" + ids[b];
        if (!phaseLinkMap.has(key))
          phaseLinkMap.set(key, { from: ids[a], to: ids[b], fragments: [] });
        phaseLinkMap.get(key).fragments.push(f.id);
      }
  }
  const phaseLinks = [...phaseLinkMap.values()].map((link) => ({
    ...link,
    support: link.fragments.length,
    linked: link.fragments.length >= opts.support,
  }));
  for (const link of phaseLinks)
    if (link.linked) parent[find(link.to)] = find(link.from);
  const blockMap = new Map();
  for (const b of bubbles) {
    const p = find(b.id);
    if (!blockMap.has(p)) blockMap.set(p, []);
    blockMap.get(p).push(b.id);
  }
  const blocks = [];
  for (const featureIDs of blockMap.values()) {
    const candidates = fragments.filter((f) =>
        featureIDs.some((id) => f.features.has(id)),
      ),
      groups = [],
      ambiguous = [],
      contradictory = [];
    const ordered = candidates
      .slice()
      .sort(
        (a, b) =>
          b.features.size - a.features.size ||
          b.guard.size - a.guard.size ||
          a.id - b.id,
      );
    for (const f of ordered) {
      if (f.conflict) {
        contradictory.push(f.id);
        continue;
      }
      const calls = new Map(
        [...f.features].filter(([id]) => featureIDs.includes(id)),
      );
      const matches = groups.filter(
        (g) =>
          [...calls.keys()].some((id) => g.calls.has(id)) &&
          compatible(g.calls, calls) &&
          compatible(g.guard, f.guard),
      );
      if (matches.length > 1) {
        ambiguous.push(f.id);
        continue;
      }
      if (matches.length === 1) {
        const g = matches[0];
        g.fragments.push(f.id);
        for (const [id, a] of calls) g.calls.set(id, a);
        for (const [p, b] of f.guard) g.guard.set(p, b);
      } else
        groups.push({
          id: groups.length,
          calls,
          guard: new Map(f.guard),
          fragments: [f.id],
        });
    }
    // Reassign ambiguous fragments only when there is one compatible group
    // after all constraints have arrived. They never act as bridging evidence.
    const shared = [];
    for (const id of ambiguous) {
      const f = fragments[id],
        matches = groups.filter(
          (g) =>
            compatible(g.calls, f.features) && compatible(g.guard, f.guard),
        );
      if (matches.length === 1) matches[0].fragments.push(id);
      else shared.push({ fragment: id, groups: matches.map((g) => g.id) });
    }
    blocks.push({
      id: blocks.length,
      features: featureIDs,
      groups: groups.map((g) => ({
        id: g.id,
        calls: [...g.calls],
        fragments: g.fragments,
        supported: g.fragments.length >= opts.support,
      })),
      shared,
      contradictory,
      total: candidates.length,
    });
  }
  const suppressedSegments = new Set(
    [...suppressedGraphNodes].map((n) => loc.get(n).unitig),
  );
  const visibleSegments = unitigs.filter((u) => !suppressedSegments.has(u.id));
  const visibleID = new Map(visibleSegments.map((u, i) => [u.id, i]));
  const structural = structuralEvidence(reads, region, opts);
  return {
    region,
    reference,
    settings: opts,
    suppressedVariants: vafFilter.hiddenVariants,
    phaseLinks,
    reads,
    fragments: fragments.map((f) => ({
      id: f.id,
      name: f.name,
      reads: f.reads,
      calls: [...f.calls],
      features: [...f.features],
      conflict: f.conflict,
    })),
    unitigs: visibleSegments.map((u, id) => ({ ...u, id })),
    links: links
      .filter((e) => visibleID.has(e.from) && visibleID.has(e.to))
      .map((e) => ({
        ...e,
        from: visibleID.get(e.from),
        to: visibleID.get(e.to),
      })),
    bubbles: bubbles.map((b) => ({
      ...b,
      callable: fragments.filter((f) => !f.conflict && f.features.has(b.id))
        .length,
      alleles: b.alleles.map((a) => ({
        id: a.id,
        seq: a.seq,
        support: a.fs.size,
        callSupport: fragments.filter(
          (f) => !f.conflict && f.features.get(b.id) === a.id,
        ).length,
        fragments: [...a.fs],
        unitigs: [
          ...new Set(
            a.nodes
              .map((n) => visibleID.get(loc.get(n).unitig))
              .filter((id) => id !== undefined),
          ),
        ],
        isReference: a.isReference,
      })),
    })),
    blocks,
    structural,
    stats: {
      records: reads.length,
      eligibleRecords: reads.filter((r) => !r.exclusion).length,
      fragments: fragments.length,
      orientedKmers: nodes.length,
      edges: edges.size,
      weakUnitigs: visibleSegments.filter((u) => u.support < opts.support)
        .length,
      complex,
      unthreaded: fragments.filter((f) => !f.features.size).length,
      issues,
    },
  };
}

export function exportGFA(result) {
  const lines = ["H\tVN:Z:1.0"];
  for (const u of result.unitigs)
    lines.push(
      `S\tu${u.id + 1}\t${u.seq}\tFC:i:${u.fragments.length}\tms:i:${u.support}`,
    );
  for (const e of result.links)
    lines.push(
      `L\tu${e.from + 1}\t${e.fromOrientation}\tu${e.to + 1}\t${e.toOrientation}\t${result.settings.k - 1}M\tFC:i:${e.support}`,
    );
  return lines.join("\n") + "\n";
}
