// Complete compatibility with every already-observed base, including reference
// observations. Missing positions add no constraint and cannot erase a clash.
export function compatibleCalls(a, b) {
  const [small, large] = a.size <= b.size ? [a, b] : [b, a];
  for (const [pos, base] of small) {
    if (base === "!") return false;
    const other = large.get(pos);
    if (other === "!" || (other !== undefined && other !== base)) return false;
  }
  return true;
}
export function mergeCalls(target, source) {
  for (const [p, b] of source) target.set(p, b);
}
export function partitionCompatible(items, getCalls) {
  const partitions = [];
  for (const item of [...items].sort(
    (a, b) => getCalls(b).size - getCalls(a).size,
  )) {
    const calls = getCalls(item);
    const internal = [...calls.values()].includes("!");
    const matches = internal
      ? []
      : partitions.filter(
          (p) => !p.internal && compatibleCalls(p.calls, calls),
        );
    // A partial read compatible with several conflicting subgroups gets its own
    // unallocated partition, rather than an arbitrary clone assignment.
    const chosen = matches.length === 1 ? matches[0] : null;
    if (chosen) {
      chosen.items.push(item);
      mergeCalls(chosen.calls, calls);
    } else
      partitions.push({
        items: [item],
        calls: new Map(calls),
        internal,
        ambiguous: matches.length > 1,
      });
  }
  return partitions;
}
