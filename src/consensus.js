import { cigarOps, normalize } from "./analysis.js";
// Summarize all members, counting overlapping mates once at each base. Minority
// errors remain in the expanded reads instead of painting a summary rainbow.
export function consensusGeometry(members, byID, reference, origin, opts) {
  const segments = [],
    votes = new Map(),
    events = new Map(),
    clips = new Map();
  const anyEligible = members.some((m) =>
    m.readIds.some((id) => !byID.get(id).excluded),
  );
  for (const m of members) {
    const calls = new Map(),
      seenEvents = new Map(),
      seenClips = new Map();
    for (const id of m.readIds) {
      const r = byID.get(id);
      let p = r.start,
        q = 0;
      const vote = !r.excluded || !anyEligible;
      for (const [n, op] of cigarOps(r.cigar)) {
        if ("M=X".includes(op)) {
          segments.push([p, p + n]);
          if (vote)
            for (let j = 0; j < n; j++) {
              const pos = p + j,
                b = r.seq[q + j],
                qual = r.qual[q + j];
              if (
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
          p += n;
          q += n;
        } else if (op === "I" || op === "D") {
          const ref =
              op === "D" ? reference.slice(p - origin, p - origin + n) : "",
            alt = op === "I" ? r.seq.slice(q, q + n) : "",
            qual =
              op === "I" ? r.qual.slice(q, q + n) : [r.qual[q - 1], r.qual[q]];
          if (
            vote &&
            p >= origin &&
            p + (op === "D" ? n : 0) <= origin + reference.length &&
            qual.length &&
            qual.every((v) => v !== undefined && v !== 255 && v >= opts.baseq)
          ) {
            const e = normalize(p, ref, alt, reference, origin);
            seenEvents.set(e.id, { ...e, op, length: n, rawPos: e.pos });
          }
          if (op === "I") q += n;
          else p += n;
        } else if (op === "N") p += n;
        else if (op === "S") {
          if (vote) {
            const c = { pos: p, seq: r.seq.slice(q, q + n) };
            seenClips.set(`${p}:${c.seq}`, c);
          }
          q += n;
        }
      }
    }
    for (const [p, alleles] of calls) {
      if (alleles.size !== 1) continue;
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
  const mismatches = [];
  for (const [pos, v] of votes) {
    const sorted = [...v].sort((a, b) => b[1] - a[1]),
      [base, count] = sorted[0],
      total = sorted.reduce((n, a) => n + a[1], 0);
    if (
      count >= opts.support &&
      count / total >= 0.7 &&
      base !== reference[pos - origin]
    )
      mismatches.push({ pos, base, support: count });
  }
  return {
    segments: merged,
    mismatches,
    events: [...events.values()].filter((e) => e.support >= opts.support),
    clips: [...clips.values()].filter((c) => c.support >= opts.support),
  };
}
