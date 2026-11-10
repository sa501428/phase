import { test } from "node:test";
import assert from "node:assert/strict";
import {
  alignmentRuns,
  prepareRecord,
  rawGeometry,
} from "../src/alignment-geometry.js";
import { phase } from "../src/analysis.js";
import { buildFamilies } from "../src/families.js";
import { buildOverview } from "../src/overview.js";
import { consensusGeometry } from "../src/consensus.js";
const opts = { mapq: 20, baseq: 20, support: 3, duplicates: false, sample: "" },
  ref = "A".repeat(100);
const read = (id, extra = {}) => ({
  id,
  name: "r" + id,
  start: 10,
  end: 20,
  cigar: "10M",
  seq: "A".repeat(10),
  qual: Array(10).fill(35),
  flags: 0,
  mq: 60,
  rg: "r",
  sample: "s",
  mate: -1,
  ...extra,
});
function analyze(reads) {
  const result = phase(reads, ref, 0, opts);
  result.families = buildFamilies(result, ref, 0, opts);
  result.overview = buildOverview(result, { start: 0, end: 100 }, opts);
  for (const group of result.overview)
    assert.doesNotThrow(() =>
      consensusGeometry(
        group.members,
        new Map(result.reads.map((r) => [r.id, r])),
        ref,
        0,
        opts,
      ),
    );
  return result;
}
test("valid omitted SEQ preserves CIGAR geometry without fabricated SNPs", () => {
  for (const seq of ["", "*"]) {
    const r = read(0, { seq, qual: [] });
    prepareRecord(r);
    assert.equal(r.displayIssue, "missing-sequence");
    assert.equal(r.sequenceUsable, false);
    assert.equal(alignmentRuns(r)[0].length, 10);
    const g = rawGeometry([r], ref, 0);
    assert.deepEqual(g.unavailableSegments, [[10, 20]]);
    assert.equal(g.mismatches.length, 0);
  }
});
test("placed unmapped records with stored sequence and empty or unavailable CIGAR do not crash", () => {
  for (const cigar of ["", "*", "2M"]) {
    const r = read(0, { flags: 4, cigar, end: 11 });
    prepareRecord(r);
    assert.deepEqual(alignmentRuns(r), []);
    assert.equal(analyze([r]).fragmentCount, 0);
  }
});
test("SEQ-omitted primary, secondary and supplementary records do not abort healthy reads", () => {
  const good = [read(0), read(1), read(2)],
    records = [
      ...good,
      read(3, { seq: "", qual: [] }),
      read(4, { seq: "", qual: [], flags: 256 }),
      read(5, { seq: "", qual: [], flags: 2048 }),
    ],
    r = analyze(records);
  assert.equal(r.fragmentCount, 3);
  assert.equal(r.recordIssues["missing-sequence"], 3);
  assert.equal(r.reads.length, 6);
  assert.equal(
    new Set(r.overview.flatMap((g) => g.members.flatMap((m) => m.readIds)))
      .size,
    6,
  );
});
test("missing-CIGAR records stay available for inspection and contribute no bases", () => {
  const r = analyze([read(0, { cigar: "", end: 11 }), read(1)]);
  assert.equal(r.fragmentCount, 1);
  assert.equal(r.recordIssues["missing-cigar"], 1);
  assert.equal(
    new Set(r.overview.flatMap((g) => g.members.flatMap((m) => m.readIds)))
      .size,
    2,
  );
});
test("genuinely inconsistent SEQ length is quarantined per record, not silently realigned", () => {
  const r = analyze([
    read(0, { seq: "AAAA", qual: Array(4).fill(35) }),
    read(1),
  ]);
  assert.equal(r.fragmentCount, 1);
  assert.equal(r.recordIssues["inconsistent-sequence-length"], 1);
  assert.equal(r.reads.length, 2);
  const g = rawGeometry(r.reads, ref, 0);
  assert.equal(g.mismatches.length, 0);
  assert.equal(g.unavailableSegments.length, 1);
});
test("hard clipping consumes no stored SEQ while soft clipping does", () => {
  const r = read(0, { cigar: "30H2S8M10H", seq: "TTAAAAAAAA", end: 18 });
  prepareRecord(r);
  assert.equal(r.displayIssue, undefined);
  assert.equal(alignmentRuns(r).find((run) => run.op === "M").queryStart, 2);
});
