// Average-linkage SNP similarity, per phase block. No ancestry direction or
// cross-block clone identity is inferred. Missing calls never count as matches.
export function snpDistance(a, b, sites) {
  let compared = 0,
    different = 0,
    sharedAlternate = 0,
    sharedReference = 0;
  for (const { index, ref } of sites) {
    const x = a[index],
      y = b[index];
    if (x === "?" || y === "?" || x === undefined || y === undefined) continue;
    compared++;
    if (x !== y) different++;
    else if (x === ref) sharedReference++;
    else sharedAlternate++;
  }
  return {
    distance: compared ? different / compared : null,
    compared,
    different,
    sharedAlternate,
    sharedReference,
  };
}
export function clusterPatterns(patterns, sites) {
  const pairs = patterns.map((a) =>
    patterns.map((b) => snpDistance(a.pattern, b.pattern, sites)),
  );
  let nodes = patterns.map((g, i) => ({
      id: "H" + g.id,
      groupID: g.id,
      leaves: [i],
      height: 0,
    })),
    serial = 0;
  while (nodes.length > 1) {
    let best = null;
    for (let i = 0; i < nodes.length; i++)
      for (let j = i + 1; j < nodes.length; j++) {
        let sum = 0,
          count = 0,
          valid = true,
          minCompared = Infinity;
        for (const a of nodes[i].leaves)
          for (const b of nodes[j].leaves) {
            const p = pairs[a][b];
            if (p.distance === null) {
              valid = false;
              break;
            }
            sum += p.distance;
            count++;
            minCompared = Math.min(minCompared, p.compared);
          }
        if (!valid) continue;
        const distance = sum / count;
        if (!best || distance < best.distance)
          best = { i, j, distance, minCompared };
      }
    if (!best) break;
    const a = nodes[best.i],
      b = nodes[best.j],
      parent = {
        id: "N" + serial++,
        left: a,
        right: b,
        leaves: [...a.leaves, ...b.leaves],
        height: best.distance,
        minCompared: best.minCompared,
      };
    nodes = nodes.filter((_, i) => i !== best.i && i !== best.j);
    nodes.push(parent);
  }
  return { roots: nodes, pairs };
}
export function cutHierarchy(roots, threshold) {
  const clusters = [];
  function visit(n) {
    if (!n.left || n.height <= threshold) clusters.push(n.leaves);
    else {
      visit(n.left);
      visit(n.right);
    }
  }
  roots.forEach(visit);
  return clusters;
}
export function buildHierarchy(result) {
  const frequencies = new Map();
  function siteCounts(index) {
    if (!frequencies.has(index)) {
      const counts = new Map();
      let total = 0;
      for (const vector of result.vectors) {
        const value = vector[index];
        if (value === "?" || value === undefined) continue;
        total++;
        counts.set(value, (counts.get(value) || 0) + 1);
      }
      frequencies.set(index, { counts, total });
    }
    return frequencies.get(index);
  }
  return result.blocks.map((block) => {
    const unique = new Map();
    block.indices.forEach((variantIndex, localIndex) => {
      const v = result.variants[variantIndex];
      if (v.ref.length === 1 && v.alt.length === 1 && !unique.has(v.pos))
        unique.set(v.pos, {
          index: localIndex,
          variantIndex,
          pos: v.pos,
          ref: v.ref,
        });
    });
    const sites = [...unique.values()];
    const informative = result.vectors.reduce(
        (n, v) => n + Number(block.indices.some((i) => v[i] !== "?")),
        0,
      ),
      assigned = block.groups.reduce((n, g) => n + g.members.length, 0);
    const groups = block.groups.map((g) => {
      const readIDs = new Set(
        g.members.flatMap((i) => result.fragments[i].reads),
      );
      const alleleFractions = [],
        seen = new Set();
      block.indices.forEach((index, i) => {
        const v = result.variants[index],
          a = g.pattern[i],
          key = v.pos + ":" + v.ref;
        if (a === "?" || a === (v.ref || "-") || seen.has(key)) return;
        seen.add(key);
        const { counts, total } = siteCounts(index),
          alt = counts.get(a) || 0;
        alleleFractions.push({
          pos: v.pos,
          ref: v.ref,
          allele: a,
          alt,
          total,
          fraction: total ? alt / total : null,
        });
      });
      return {
        id: g.id,
        pattern: g.pattern,
        fragments: g.members.length,
        records: readIDs.size,
        informative,
        assignedShare: informative ? g.members.length / informative : null,
        ambiguousCompatible: g.compatible,
        alleleFractions,
      };
    });
    return {
      blockID: block.id,
      start: block.start,
      end: block.end,
      groups,
      sites,
      informative,
      assigned,
      ambiguous: block.ambiguous,
      conflicting: block.conflicting,
      ...clusterPatterns(groups, sites),
    };
  });
}
