import { test } from "node:test";
import assert from "node:assert/strict";
import { alignmentRuns, rawGeometry } from "../src/alignment-geometry.js";
import { consensusGeometry } from "../src/consensus.js";
const opts = { baseq: 20, support: 3 };
function consensus(reads, ref, origin) {
  return consensusGeometry(
    reads.map((r) => ({ readIds: [r.id] })),
    new Map(reads.map((r) => [r.id, r])),
    ref,
    origin,
    opts,
  );
}
test("collapsed insertion and deletion glyphs keep the same BAM anchors as expanded reads in a repeat", () => {
  const ref = "A".repeat(80),
    seq = "A".repeat(8) + "C" + "A".repeat(4) + "G" + "AA";
  const reads = Array.from({ length: 3 }, (_, id) => ({
    id,
    start: 103,
    end: 119,
    cigar: "5M1I5M1D5M",
    seq,
    qual: Array(16).fill(35),
    flags: id % 2 ? 16 : 0,
  }));
  const raw = rawGeometry(reads, ref, 100),
    collapsed = consensus(reads, ref, 100);
  assert.deepEqual(
    raw.events.map((e) => [e.op, e.rawPos]),
    [
      ["I", 108],
      ["D", 113],
    ],
  );
  assert.deepEqual(
    collapsed.events.map((e) => [e.op, e.rawPos]),
    [
      ["I", 108],
      ["D", 113],
    ],
  );
  assert.deepEqual(
    collapsed.events.map((e) => e.pos),
    [100, 100],
  ); // analysis-normalized positions intentionally differ
  assert.deepEqual(
    collapsed.mismatches.map((m) => [m.pos, m.base]),
    [
      [110, "C"],
      [116, "G"],
    ],
  );
  assert.deepEqual(
    raw.mismatches,
    collapsed.mismatches.map(({ support, ...m }) => m),
  );
});
test("SNPs after soft clips, insertions, deletions and skips use identical genomic columns", () => {
  const ref = "A".repeat(40),
    seq = "TTT" + "AAAA" + "TT" + "ACA" + "AA" + "AGA" + "T";
  const reads = Array.from({ length: 3 }, (_, id) => ({
    id,
    start: 100,
    end: 119,
    cigar: "2H3S4M2I3M2D2M5N3M1P1S",
    seq,
    qual: Array(seq.length).fill(35),
    flags: id % 2 ? 16 : 0,
  }));
  const runs = alignmentRuns(reads[0]);
  assert.deepEqual(
    runs
      .filter((r) => "M=X".includes(r.op))
      .map((r) => [r.refStart, r.queryStart]),
    [
      [100, 3],
      [104, 9],
      [109, 12],
      [116, 14],
    ],
  );
  const raw = rawGeometry(reads, ref, 100),
    collapsed = consensus(reads, ref, 100);
  assert.deepEqual(
    raw.mismatches.map((m) => [m.pos, m.base]),
    [
      [105, "C"],
      [117, "G"],
    ],
  );
  assert.deepEqual(
    collapsed.mismatches.map((m) => [m.pos, m.base]),
    [
      [105, "C"],
      [117, "G"],
    ],
  );
});
test("repeat-equivalent raw anchors are not moved or pooled into a fictitious display column", () => {
  const ref = "A".repeat(30),
    reads = Array.from({ length: 3 }, (_, id) => ({
      id,
      start: 100,
      end: 110,
      cigar: `${3 + id}M1I${7 - id}M`,
      seq: "A".repeat(11),
      qual: Array(11).fill(35),
    }));
  const raw = rawGeometry(reads, ref, 100),
    collapsed = consensus(reads, ref, 100);
  assert.deepEqual(
    raw.events.map((e) => e.rawPos),
    [103, 104, 105],
  );
  assert.equal(collapsed.events.length, 0); // none of the raw placements has three independent votes
});
test("malformed display CIGARs fail instead of shifting subsequent SNP columns", () => {
  assert.throws(
    () => alignmentRuns({ start: 0, cigar: "2M1I2M", seq: "AAAA" }),
    /Sequence length/,
  );
});
