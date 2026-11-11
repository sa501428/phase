import { hiddenBase, hiddenEvent } from "./vaf-filter.js";
import { cigarOps } from "./cigar.js";
export const hasSequence = (read) =>
  read.sequenceUsable !== false &&
  typeof read.seq === "string" &&
  read.seq !== "" &&
  read.seq !== "*";
export function prepareRecord(read) {
  read.alignmentUsable = true;
  read.sequenceUsable = true;
  delete read.displayIssue;
  if (read.flags & 4) {
    read.alignmentUsable = false;
    read.sequenceUsable = false;
    read.displayIssue = "unmapped";
    return read;
  }
  if (!read.cigar || read.cigar === "*") {
    read.alignmentUsable = false;
    read.sequenceUsable = false;
    read.displayIssue = "missing-cigar";
    return read;
  }
  let ops;
  try {
    ops = cigarOps(read.cigar);
  } catch {
    read.alignmentUsable = false;
    read.sequenceUsable = false;
    read.displayIssue = "invalid-cigar";
    return read;
  }
  if (!read.seq || read.seq === "*") {
    read.sequenceUsable = false;
    read.displayIssue = "missing-sequence";
    return read;
  }
  const length = ops.reduce(
    (n, [size, op]) => n + ("M=XIS".includes(op) ? size : 0),
    0,
  );
  if (length !== read.seq.length) {
    read.sequenceUsable = false;
    read.displayIssue = "inconsistent-sequence-length";
  }
  return read;
}
// BAM SEQ is already oriented for its alignment. All pileup rows use the same
// zero-based genomic axis, regardless of clipping, strand, or consensus state.
export function alignmentRuns(read) {
  if (
    read.flags & 4 ||
    read.alignmentUsable === false ||
    !read.cigar ||
    read.cigar === "*"
  )
    return [];
  let refStart = read.start,
    queryStart = 0;
  const runs = [];
  for (const [length, op] of cigarOps(read.cigar)) {
    runs.push({ op, length, refStart, queryStart });
    if ("M=XDN".includes(op)) refStart += length;
    if ("M=XIS".includes(op)) queryStart += length;
  }
  if (hasSequence(read) && queryStart !== read.seq.length)
    throw Error("Sequence length does not match display CIGAR");
  return runs;
}
export function rawGeometry(reads, reference, origin) {
  const segments = [],
    mismatches = new Map(),
    events = new Map(),
    clips = new Map();
  const unavailableSegments = [];
  for (const r of reads)
    for (const { op, length: n, refStart: p, queryStart: q } of alignmentRuns(
      r,
    )) {
      if ("M=X".includes(op)) {
        if (!hasSequence(r)) {
          unavailableSegments.push([p, p + n]);
          continue;
        }
        segments.push([p, p + n]);
        for (let j = 0; j < n; j++) {
          const pos = p + j,
            base = r.seq[q + j];
          if (
            !hiddenBase(r, q + j) &&
            pos >= origin &&
            pos < origin + reference.length &&
            base &&
            base !== reference[pos - origin]
          )
            mismatches.set(`${pos}:${base}`, { pos, base });
        }
      } else if (
        (op === "I" || op === "D" || op === "N") &&
        !hiddenEvent(r, op, p)
      )
        events.set(`${p}:${op}:${n}`, { pos: p, rawPos: p, op, length: n });
      else if (op === "S" && hasSequence(r)) {
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
    unavailableSegments,
    mismatches: [...mismatches.values()],
    events: [...events.values()],
    clips: [...clips.values()],
  };
}
