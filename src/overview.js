// A compact view of inferred LOCAL groups. Partial/low-quality observations do
// not create display rows merely because their endpoints or error bases differ.
export function buildOverview(result, region, opts) {
  const byID = new Map(result.reads.map((r) => [r.id, r])),
    vectorsByRead = new Map();
  for (const f of result.fragments)
    for (const id of f.reads) vectorsByRead.set(id, result.vectors[f.id]);
  const templates = result.families.flatMap((f) =>
    f.members.map((m) => ({
      ...m,
      displayVector: f.observed?.map((a, i) =>
        a === "ALT"
          ? result.variants[i].alt || "-"
          : a === "REF"
            ? result.variants[i].ref || "-"
            : a,
      ),
    })),
  );
  const rows = [],
    supportedSV = new Set(
      (result.structural?.events || [])
        .filter(
          (e) =>
            e.eligibleFragments >= opts.support && e.kind !== "ambiguous-split",
        )
        .map((e) => e.id),
    );
  const scopes = result.blocks.map((b, i, blocks) => ({
    block: b,
    start: i
      ? blocks[i - 1].end <= b.start
        ? Math.floor((blocks[i - 1].end + b.start) / 2)
        : Math.min(blocks[i - 1].start, b.start)
      : region.start,
    end:
      i + 1 < blocks.length
        ? b.end <= blocks[i + 1].start
          ? Math.floor((b.end + blocks[i + 1].start) / 2)
          : Math.max(b.end, blocks[i + 1].end)
        : region.end,
  }));
  if (!scopes.length) {
    const components = [];
    for (const m of [...templates].sort((a, b) => a.start - b.start)) {
      const start = Math.max(region.start, m.start),
        end = Math.min(region.end, m.end);
      if (end <= start) continue;
      const prev = components.at(-1);
      if (prev && start < prev.end) prev.end = Math.max(prev.end, end);
      else components.push({ start, end });
    }
    scopes.push(...components);
  }
  for (const scope of scopes) {
    const buckets = new Map(),
      block = scope.block;
    function bucket(key, label, category, pattern, m) {
      if (!buckets.has(key))
        buckets.set(key, {
          id: rows.length + buckets.size,
          kind: "overview",
          label,
          category,
          groups: [label],
          scope: { start: scope.start, end: scope.end },
          blockID: block?.id,
          pattern,
          indices: block?.indices || [],
          members: [],
          start: scope.start,
          end: scope.end,
        });
      buckets.get(key).members.push(m);
    }
    for (const m of templates) {
      if (!m.segments.some(([s, e]) => s < scope.end && e > scope.start))
        continue;
      const rr = m.readIds.map((id) => byID.get(id)),
        eligible = rr.some((r) => !r.excluded),
        eligibleVector = m.readIds
          .map((id) => vectorsByRead.get(id))
          .find(Boolean);
      const duplicateOnly =
        !eligible &&
        rr.every(
          (r) =>
            !(r.flags & (4 | 256 | 512)) && r.mq !== 255 && r.mq >= opts.mapq,
        ) &&
        rr.some((r) => r.flags & 1024);
      if (!eligible && !duplicateOnly) {
        bucket("excluded", "Excluded evidence", "excluded", null, m);
        continue;
      }
      if (!block) {
        const sv = [
          ...new Set(
            rr
              .flatMap((r) => r.structuralLinks || [])
              .filter((id) => supportedSV.has(id)),
          ),
        ].sort((a, b) => a - b);
        bucket(
          "coverage:" + sv.join(","),
          sv.length
            ? "Coverage · SV " + sv.join(",")
            : supportedSV.size
              ? "Coverage · SV linkage unknown"
              : "Coverage consensus",
          "coverage",
          null,
          m,
        );
        continue;
      }
      const vector =
          eligibleVector || m.displayVector || result.variants.map(() => "?"),
        a = block.indices.map((i) => vector[i] ?? "?"),
        informative = a.some((v) => v !== "?");
      if (!informative) {
        bucket(
          "uninformative",
          `B${block.id} · uninformative`,
          "uninformative",
          null,
          m,
        );
        continue;
      }
      const matches = block.groups.filter((g) =>
        a.every(
          (v, i) => v === "?" || g.pattern[i] === "?" || v === g.pattern[i],
        ),
      );
      const shared =
        matches.length === 1 &&
        a.some((v, i) => v !== "?" && matches[0].pattern[i] !== "?");
      if (matches.length === 1 && shared) {
        const g = matches[0],
          sv = [
            ...new Set(
              rr
                .flatMap((r) => r.structuralLinks || [])
                .filter((id) => supportedSV.has(id)),
            ),
          ].sort((a, b) => a - b),
          suffix = sv.length ? " · SV " + sv.join(",") : "";
        bucket(
          "H" + g.id + ":" + sv.join(","),
          `B${block.id} H${g.id}${suffix}`,
          "haplotype",
          g.pattern,
          m,
        );
      } else if (matches.length > 1)
        bucket("ambiguous", `B${block.id} · ambiguous`, "ambiguous", null, m);
      else
        bucket(
          "unresolved",
          `B${block.id} · unresolved`,
          "unresolved",
          null,
          m,
        );
    }
    for (const f of buckets.values()) {
      f.members.sort((a, b) => a.start - b.start || a.id - b.id);
      f.count = f.members.length;
      f.recordCount = f.members.reduce((n, m) => n + m.readIds.length, 0);
      f.excluded = f.members.filter((m) =>
        m.readIds.every((id) => byID.get(id).excluded),
      ).length;
      f.duplicates = f.members.filter((m) =>
        m.readIds.some((id) => byID.get(id).flags & 1024),
      ).length;
      f.structuralLinks = [
        ...new Set(
          f.members.flatMap((m) =>
            m.readIds.flatMap((id) => byID.get(id).structuralLinks || []),
          ),
        ),
      ];
      rows.push(f);
    }
  }
  return rows;
}
