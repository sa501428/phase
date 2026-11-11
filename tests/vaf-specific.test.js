import { test } from "node:test";
import assert from "node:assert/strict";
import { analyzeReads } from "../specific/helpers.js";
const variants = [
  {
    id: 0,
    chr: "chr1",
    pos: 30,
    ref: "A",
    alt: "C",
    refLength: 1,
    label: "rare",
  },
  {
    id: 1,
    chr: "chr1",
    pos: 40,
    ref: "A",
    alt: "G",
    refLength: 1,
    label: "common",
  },
];
const opts = { mapq: 20, baseq: 20, duplicates: false, helpers: false };
function reads(n) {
  return Array.from({ length: n }, (_, id) => {
    let seq = Array(100).fill("A");
    if (!id) seq[30] = "C";
    if (id < 10) seq[40] = "G";
    return {
      id,
      name: "r" + id,
      chr: "chr1",
      start: 0,
      seq: seq.join(""),
      cigar: "100M",
      qual: Array(100).fill(35),
      mq: 60,
      flags: 0,
    };
  });
}
test("specific viewer removes rare requested variants from phase pairs and reindexes remaining evidence", () => {
  const r = analyzeReads(reads(51), variants, opts);
  assert.equal(r.requestedCount, 1);
  assert.equal(r.variants[0].label, "common");
  assert.equal(r.variants[0].id, 0);
  assert.equal(r.sites[0].alt, 10);
  assert.equal(r.pairs.length, 0);
  assert.equal(r.suppressedVariants[0].label, "rare");
  assert.ok([...r.collector.fragments.values()].every((f) => !f.calls.has(1)));
  const all = analyzeReads(reads(51), variants, { ...opts, minVaf: 0 });
  assert.equal(all.requestedCount, 2);
  assert.equal(all.pairs[0].both, 1);
  assert.equal(analyzeReads(reads(50), variants, opts).requestedCount, 2);
});
test("specific viewer handles all variants suppressed and rejects invalid cutoffs", () => {
  const r = analyzeReads(reads(51), variants, { ...opts, minVaf: 0.9 });
  assert.equal(r.requestedCount, 0);
  assert.equal(r.variants.length, 0);
  assert.equal(r.pairs.length, 0);
  for (const v of [-0.1, 1.1, NaN])
    assert.throws(
      () => analyzeReads(reads(2), variants, { ...opts, minVaf: v }),
      /VAF/,
    );
});
