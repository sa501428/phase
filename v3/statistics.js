import { clusterPatterns, cutHierarchy } from "../src/hierarchy.js";
export function siteCounts(model, ids, pos, insertion = false) {
  const counts = new Map();
  let total = 0,
    unknown = 0,
    conflicting = 0;
  for (const id of model.activeIDs(ids)) {
    const e = model.evidence(id),
      values = (insertion ? e.insertions : e.bases).get(pos);
    if (!values) {
      unknown++;
      continue;
    }
    if (values.size !== 1) {
      conflicting++;
      continue;
    }
    const value = [...values][0];
    counts.set(value, (counts.get(value) || 0) + 1);
    total++;
  }
  return { counts, total, unknown, conflicting };
}
export function groupStatistics(model, g) {
  const summary = model.groupSummary(g),
    all = model.fragments.map((f) => f.id),
    total = model.activeIDs(all).length,
    members = model.activeIDs(g.members),
    variants = [];
  for (const c of summary.bases) {
    if (c.conflict) continue;
    const alt = c.values[0].value,
      ref = model.reference[c.pos - model.region.start];
    if (!ref || alt === ref) continue;
    const global = siteCounts(model, all, c.pos),
      local = siteCounts(model, g.members, c.pos);
    variants.push({
      pos: c.pos,
      ref,
      alt,
      kind: alt === "-" ? "deletion" : "SNV",
      altCount: global.counts.get(alt) || 0,
      total: global.total,
      localAlt: local.counts.get(alt) || 0,
      localTotal: local.total,
    });
  }
  for (const c of summary.insertions) {
    if (c.conflict) continue;
    const alt = c.values[0].value,
      global = siteCounts(model, all, c.pos, true),
      local = siteCounts(model, g.members, c.pos, true);
    variants.push({
      pos: c.pos,
      ref: "",
      alt,
      kind: "insertion",
      altCount: global.counts.get(alt) || 0,
      total: global.total,
      localAlt: local.counts.get(alt) || 0,
      localTotal: local.total,
    });
  }
  return {
    members: members.length,
    assigned: g.members.length,
    hidden: g.members.filter((id) => model.hidden.has(id)).length,
    total,
    share: total ? members.length / total : null,
    variants,
    summary,
  };
}
export function manualHierarchy(model, threshold = 0.25) {
  if (model.groups.length > 80)
    return {
      error:
        "Hierarchy supports up to 80 constructed haplotypes. All assignments and reads are retained.",
    };
  const calls = model.groups.map(
      (g) =>
        new Map(
          model
            .groupSummary(g)
            .bases.filter(
              (c) => !c.conflict && /^[ACGT]$/.test(c.values[0].value),
            )
            .map((c) => [c.pos, c.values[0].value]),
        ),
    ),
    candidate = new Map();
  for (const m of calls)
    for (const [p, b] of m) {
      if (!candidate.has(p)) candidate.set(p, new Set());
      candidate.get(p).add(b);
    }
  const positions = [...candidate]
      .filter(
        ([p, bs]) =>
          bs.size > 1 ||
          [...bs].some((b) => b !== model.reference[p - model.region.start]),
      )
      .map(([p]) => p)
      .sort((a, b) => a - b),
    sites = positions.map((p, index) => ({
      index,
      ref: model.reference[p - model.region.start],
      pos: p,
    })),
    groups = model.groups.map((g, i) => ({
      id: g.id,
      name: g.name,
      pattern: positions.map((p) => calls[i].get(p) ?? "?"),
    })),
    tree = clusterPatterns(groups, sites);
  return {
    ...tree,
    groups,
    sites,
    clusters: cutHierarchy(tree.roots, threshold),
  };
}
