import { cigarOps, normalize } from "./analysis.js";

// Display families never change analysis counts. Unknown/discordant evidence
// stays separate, and fragment endpoints do not define an allele family.
export function buildFamilies(result, reference, origin, opts) {
  const templates = new Map(),
    assignment = new Map();
  for (const f of result.fragments)
    for (const id of f.reads)
      assignment.set(id, f.groups.join(", ") || "Unassigned");
  for (const r of result.reads) {
    const key = JSON.stringify([
      r.sample,
      r.rg,
      r.name,
      r.flags & (1 | 2048) ? "template" : r.id,
    ]);
    if (!templates.has(key)) templates.set(key, []);
    templates.get(key).push(r);
  }
  const components = [];
  for (const rr of [...templates.values()].sort(
    (a, b) =>
      Math.min(...a.map((r) => r.start)) - Math.min(...b.map((r) => r.start)),
  )) {
    const start = Math.min(...rr.map((r) => r.start)),
      end = Math.max(...rr.map((r) => r.end)),
      last = components.at(-1);
    if (last && start < last[1]) last[1] = Math.max(last[1], end);
    else components.push([start, end]);
  }
  function componentAt(position) {
    let lo = 0,
      hi = components.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (components[mid][0] <= position) lo = mid + 1;
      else hi = mid;
    }
    return lo - 1;
  }
  const families = new Map();
  for (const reads of templates.values()) {
    const differences = new Map(),
      calls = new Map(),
      segments = [],
      events = [],
      clips = [];
    let uncertain = false;
    let conflictingSNP = false;
    const rawCalls = new Map();
    const recordEvidence = [];
    for (const r of reads) {
      if (r.start < origin || r.end > origin + reference.length)
        uncertain = true;
      const rSegments = [],
        rEvents = [];
      let p = r.start,
        q = 0;
      for (const [n, op] of cigarOps(r.cigar)) {
        if ("M=X".includes(op)) {
          segments.push([p, p + n]);
          rSegments.push([p, p + n]);
          for (let i = 0; i < n; i++) {
            const pos = p + i,
              b = r.seq[q + i],
              ref = reference[pos - origin],
              quality = r.qual[q + i];
            if (!ref) continue;
            if ("ACGT".includes(b)) {
              if (rawCalls.has(pos) && rawCalls.get(pos) !== b)
                conflictingSNP = true;
              rawCalls.set(pos, b);
            }
            const qualified =
              quality !== undefined && quality !== 255 && quality >= opts.baseq;
            const prev = calls.get(pos);
            if (prev && prev !== b) uncertain = true;
            if (qualified) calls.set(pos, b);
            if (b !== ref) {
              differences.set(`${pos}:S:${b}`, { pos, base: b });
              if (!qualified) uncertain = true;
            }
          }
          p += n;
          q += n;
        } else if (op === "I" || op === "D") {
          const ref =
              op === "D" ? reference.slice(p - origin, p - origin + n) : "",
            alt = op === "I" ? r.seq.slice(q, q + n) : "";
          const e = normalize(p, ref, alt, reference, origin);
          events.push({ ...e, op, length: n, rawPos: p });
          rEvents.push(e);
          if (op === "I") q += n;
          else p += n;
        } else if (op === "S") {
          clips.push({
            pos: p,
            seq: r.seq.slice(q, q + n),
            side: q === 0 ? "left" : "right",
          });
          q += n;
        } else if (op === "N") {
          events.push({ pos: p, op, length: n, id: `${p}:N:${n}` });
          p += n;
        } else if (op === "H") {
          uncertain = true;
        }
      }
      recordEvidence.push({ segments: rSegments, events: rEvents });
    }
    for (const e of events) {
      if (e.op === "N") continue;
      if (
        recordEvidence.some(
          (r) =>
            !r.events.some((v) => v.id === e.id) &&
            r.segments.some(
              ([s, end]) =>
                s < e.pos && end > e.pos + Math.max(1, e.ref.length),
            ),
        )
      )
        uncertain = true;
    }
    const snps = [...differences.entries()].sort(([a], [b]) =>
      a.localeCompare(b),
    );
    const eventIDs = [...new Set(events.map((e) => e.id))].sort();
    const clipIDs = [...new Set(clips.map((c) => JSON.stringify(c)))].sort();
    const observed = result.variants.map((v) => {
      const alt = events.some((e) => e.id === v.id);
      if (alt) return "ALT";
      if (!v.ref.length || !v.alt.length)
        return segments.some(
          ([s, e]) => s < v.pos && e > v.pos + Math.max(1, v.ref.length),
        ) &&
          calls.has(v.pos - 1) &&
          calls.has(v.pos + Math.max(1, v.ref.length))
          ? "REF"
          : "?";
      return Array.from(
        { length: v.ref.length },
        (_, i) => calls.get(v.pos + i) || "?",
      ).join("");
    });
    // Reference-only fragments in nonoverlapping regions must not form a family.
    const blockScope = result.blocks
      .filter((b) => segments.some(([s, e]) => s < b.end && e > b.start))
      .map((b) => b.id);
    const key = JSON.stringify([
      reads[0].sample,
      reads[0].rg,
      snps.map(([k]) => k),
      eventIDs,
      clipIDs,
      [...new Set(reads.flatMap((r) => r.structuralLinks || []))].sort(
        (a, b) => a - b,
      ),
      observed,
      blockScope.length
        ? blockScope
        : ["coverage", componentAt(Math.min(...reads.map((r) => r.start)))],
      conflictingSNP ? reads.map((r) => r.id) : null,
      uncertain
        ? reads
            .map((r) => [
              r.start,
              r.end,
              r.cigar,
              r.seq,
              r.qual,
              r.flags & ~1024,
            ])
            .sort((a, b) => JSON.stringify(a).localeCompare(JSON.stringify(b)))
        : null,
    ]);
    let family = families.get(key);
    if (!family) {
      family = {
        id: families.size,
        observed,
        discordantPair: conflictingSNP,
        members: [],
        segments: [],
        mismatches: snps.map(([, v]) => v),
        events: [...new Map(events.map((e) => [e.id, e])).values()],
        clips,
        groups: new Set(),
        start: Infinity,
        end: -Infinity,
        excluded: 0,
        duplicates: 0,
      };
      families.set(key, family);
    }
    const member = {
      id: reads[0].id,
      readIds: reads.map((r) => r.id),
      start: Math.min(...reads.map((r) => r.start)),
      end: Math.max(...reads.map((r) => r.end)),
      paired: reads.some((r) => r.flags & 1),
      matePresent:
        reads.some((r) => r.flags & 64) && reads.some((r) => r.flags & 128),
      segments,
      events,
      clips,
      mismatches: snps.map(([, v]) => v),
    };
    family.members.push(member);
    family.segments.push(...segments);
    family.start = Math.min(family.start, member.start);
    family.end = Math.max(family.end, member.end);
    if (reads.every((r) => r.excluded)) family.excluded++;
    if (reads.some((r) => r.flags & 1024)) family.duplicates++;
    for (const r of reads)
      family.groups.add(assignment.get(r.id) || "Excluded");
  }
  return [...families.values()].map((f) => {
    const merged = [];
    for (const [s, e] of f.segments.sort((a, b) => a[0] - b[0])) {
      const last = merged.at(-1);
      if (last && s <= last[1]) last[1] = Math.max(last[1], e);
      else merged.push([s, e]);
    }
    return {
      ...f,
      segments: merged,
      groups: [...f.groups].sort(),
      count: f.members.length,
      recordCount: f.members.reduce((n, m) => n + m.readIds.length, 0),
      structuralLinks: [
        ...new Set(
          f.members.flatMap((m) =>
            m.readIds.flatMap((id) => result.reads[id]?.structuralLinks || []),
          ),
        ),
      ],
    };
  });
}
function compareGroupRows(a, b) {
  if (a.family.kind === "overview" && b.family.kind === "overview") {
    const rank = {
      haplotype: 0,
      coverage: 0,
      ambiguous: 1,
      uninformative: 2,
      unresolved: 3,
      excluded: 4,
    };
    return (
      (a.family.blockID || 0) - (b.family.blockID || 0) ||
      a.family.scope.start - b.family.scope.start ||
      rank[a.family.category] - rank[b.family.category] ||
      a.group.localeCompare(b.group)
    );
  }
  return a.group.localeCompare(b.group) || a.start - b.start;
}
export function layoutFamilies(
  families,
  reads,
  {
    start,
    end,
    layer = "all",
    collapsed = true,
    expanded = new Set(),
    sort = "group",
  },
) {
  const byID = new Map(reads.map((r) => [r.id, r])),
    out = [];
  for (const f of families) {
    const members = f.members
      .filter((m) => m.start < end && m.end > start)
      .map((m) => ({
        ...m,
        readIds: m.readIds.filter(
          (id) =>
            layer === "all" ||
            (layer === "excluded"
              ? byID.get(id).excluded
              : !byID.get(id).excluded),
        ),
      }))
      .filter((m) => m.readIds.length);
    if (!members.length) continue;
    // Membership is never sampled; viewport rendering only determines painted rows.
    const group = f.groups.join(", "),
      base = { family: f, group, start: f.start, end: f.end };
    if (collapsed && (members.length > 1 || f.kind === "overview")) {
      out.push({ ...base, type: "family", members, count: members.length });
      if (!expanded.has(f.id)) continue;
    }
    for (const m of members)
      out.push({
        ...base,
        type: "fragment",
        member: m,
        start: m.start,
        end: m.end,
        reads: m.readIds.map((id) => byID.get(id)),
        count: 1,
      });
  }
  // Sort complete families together, preserving their header/child order.
  const buckets = new Map();
  for (const row of out) {
    if (!buckets.has(row.family.id)) buckets.set(row.family.id, []);
    buckets.get(row.family.id).push(row);
  }
  return [...buckets.values()]
    .sort((a, b) =>
      sort === "group"
        ? compareGroupRows(a[0], b[0])
        : sort === "quality"
          ? Math.max(
              ...b[0].family.members.flatMap((m) =>
                m.readIds.map((id) => byID.get(id).mq),
              ),
            ) -
            Math.max(
              ...a[0].family.members.flatMap((m) =>
                m.readIds.map((id) => byID.get(id).mq),
              ),
            )
          : sort === "strand"
            ? (byID.get(a[0].family.members[0].readIds[0]).flags & 16) -
                (byID.get(b[0].family.members[0].readIds[0]).flags & 16) ||
              a[0].start - b[0].start
            : a[0].start - b[0].start,
    )
    .flat();
}
