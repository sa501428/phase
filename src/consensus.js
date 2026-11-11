import { hiddenBase, hiddenEvent } from "./vaf-filter.js";
import { normalize } from "./analysis.js";
import { hasSequence, alignmentRuns } from "./alignment-geometry.js";
// Summarize all members, counting overlapping mates once at each base.
// Contradictory calls are flagged rather than replaced by a majority allele.
export function consensusGeometry(members, byID, reference, origin, opts) {
  const segments = [],
    votes = new Map(),
    events = new Map(),
    clips = new Map();
  const unavailableSegments = [];
  const withinFragmentConflicts = new Map();
  const anyEligible = members.some((m) =>
    m.readIds.some((id) => !byID.get(id).excluded),
  );
  for (const m of members) {
    const calls = new Map(),
      seenEvents = new Map(),
      seenClips = new Map();
    for (const id of m.readIds) {
      const r = byID.get(id);
      const vote = hasSequence(r) && (!r.excluded || !anyEligible);
      for (const { op, length: n, refStart: p, queryStart: q } of alignmentRuns(
        r,
      )) {
        if ("M=X".includes(op)) {
          if (!hasSequence(r)) unavailableSegments.push([p, p + n]);
          else segments.push([p, p + n]);
          if (vote)
            for (let j = 0; j < n; j++) {
              const pos = p + j,
                b = r.seq[q + j],
                qual = r.qual[q + j];
              if (
                hiddenBase(r, q + j) ||
                pos < origin ||
                pos >= origin + reference.length ||
                qual === undefined ||
                qual === 255 ||
                qual < opts.baseq ||
                !"ACGT".includes(b)
              )
                continue;
              if (!calls.has(pos)) calls.set(pos, new Set());
              calls.get(pos).add(b);
            }
        } else if (op === "I" || op === "D") {
          const ref =
              op === "D" ? reference.slice(p - origin, p - origin + n) : "",
            alt = op === "I" ? r.seq.slice(q, q + n) : "",
            qual =
              op === "I" ? r.qual.slice(q, q + n) : [r.qual[q - 1], r.qual[q]];
          if (
            !hiddenEvent(r, op, p) &&
            vote &&
            p >= origin &&
            p + (op === "D" ? n : 0) <= origin + reference.length &&
            qual.length &&
            qual.every((v) => v !== undefined && v !== 255 && v >= opts.baseq)
          ) {
            const e = normalize(p, ref, alt, reference, origin);
            // Normalization identifies equivalent alleles for analysis; a
            // pileup glyph must stay at its observed BAM CIGAR anchor.
            seenEvents.set(`${e.id}@${p}`, { ...e, op, length: n, rawPos: p });
          }
        } else if (op === "S") {
          if (vote) {
            const c = { pos: p, seq: r.seq.slice(q, q + n) };
            seenClips.set(`${p}:${c.seq}`, c);
          }
        }
      }
    }
    for (const [p, alleles] of calls) {
      if (alleles.size !== 1) {
        withinFragmentConflicts.set(p, { pos: p, alleles: [...alleles] });
        continue;
      }
      const b = [...alleles][0];
      if (!votes.has(p)) votes.set(p, new Map());
      const v = votes.get(p);
      v.set(b, (v.get(b) || 0) + 1);
    }
    for (const [id, e] of seenEvents) {
      if (!events.has(id)) events.set(id, { ...e, support: 0 });
      events.get(id).support++;
    }
    for (const [id, c] of seenClips) {
      if (!clips.has(id)) clips.set(id, { ...c, support: 0 });
      clips.get(id).support++;
    }
  }
  const merged = [];
  for (const [s, e] of segments.sort((a, b) => a[0] - b[0])) {
    const last = merged.at(-1);
    if (last && s <= last[1]) last[1] = Math.max(last[1], e);
    else merged.push([s, e]);
  }
  const mismatches = [],
    conflicts = [];
  for (const [pos, v] of votes) {
    if (v.size > 1) {
      conflicts.push({ pos, alleles: [...v.keys()] });
      continue;
    }
    const [base, count] = [...v][0];
    if (count >= 1 && base !== reference[pos - origin])
      mismatches.push({ pos, base, support: count });
  }
  return {
    segments: merged,
    unavailableSegments,
    conflicts: [
      ...new Map(
        [...conflicts, ...withinFragmentConflicts.values()].map((c) => [
          c.pos,
          c,
        ]),
      ).values(),
    ],
    mismatches,
    events: [...events.values()].filter((e) => e.support >= opts.support),
    clips: [...clips.values()].filter((c) => c.support >= opts.support),
  };
}
