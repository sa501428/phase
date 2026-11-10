import { cigarOps } from "./analysis.js";
// BAM SEQ is already oriented for its alignment. All pileup rows use the same
// zero-based genomic axis, regardless of clipping, strand, or consensus state.
export function alignmentRuns(read) {
  let refStart = read.start,
    queryStart = 0;
  const runs = [];
  for (const [length, op] of cigarOps(read.cigar)) {
    runs.push({ op, length, refStart, queryStart });
    if ("M=XDN".includes(op)) refStart += length;
    if ("M=XIS".includes(op)) queryStart += length;
  }
  if (queryStart !== read.seq.length)
    throw Error("Sequence length does not match display CIGAR");
  return runs;
}
export function rawGeometry(reads, reference, origin) {
  const segments = [],
    mismatches = new Map(),
    events = new Map(),
    clips = new Map();
  for (const r of reads)
    for (const { op, length: n, refStart: p, queryStart: q } of alignmentRuns(
      r,
    )) {
      if ("M=X".includes(op)) {
        segments.push([p, p + n]);
        for (let j = 0; j < n; j++) {
          const pos = p + j,
            base = r.seq[q + j];
          if (
            pos >= origin &&
            pos < origin + reference.length &&
            base !== reference[pos - origin]
          )
            mismatches.set(`${pos}:${base}`, { pos, base });
        }
      } else if (op === "I" || op === "D" || op === "N")
        events.set(`${p}:${op}:${n}`, { pos: p, rawPos: p, op, length: n });
      else if (op === "S") {
        const clip = { pos: p, seq: r.seq.slice(q, q + n) };
        clips.set(`${p}:${clip.seq}`, clip);
      }
    }
  const merged = [];
  for (const [s, e] of segments.sort((a, b) => a[0] - b[0])) {
    const last = merged.at(-1);
    if (last && s <= last[1]) last[1] = Math.max(last[1], e);
    else merged.push([s, e]);
  }
  return {
    segments: merged,
    mismatches: [...mismatches.values()],
    events: [...events.values()],
    clips: [...clips.values()],
  };
}
