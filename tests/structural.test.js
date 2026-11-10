import { test } from "node:test";
import assert from "node:assert/strict";
import { parseSA, structuralEvidence } from "../src/structural.js";
import { phase } from "../src/analysis.js";
import { buildFamilies } from "../src/families.js";
const region = { chr: "chr1", start: 1000, end: 2000 },
  opts = { mapq: 20, baseq: 20, support: 3, duplicates: false, sample: "" };
const read = (id, extra = {}) => ({
  id,
  name: "r" + id,
  rg: "a",
  sample: "s",
  chr: "chr1",
  start: 1100,
  end: 1150,
  cigar: "50M50S",
  seq: "A".repeat(100),
  qual: Array(100).fill(35),
  flags: 0,
  mq: 60,
  mate: -1,
  ...extra,
});
test("SA coordinates are converted once and malformed entries remain nonfatal", () => {
  const p = parseSA("chr2,5001,+,50S50M,60,0;chr3,bad,-,50M,60,0;");
  assert.equal(p.alignments[0].pos, 5000);
  assert.equal(p.alignments[0].qStart, 50);
  assert.equal(p.invalid, 1);
  assert.equal(parseSA("x".repeat(65537)).invalid, 1);
});
test("cross-chromosome complementary split alignment is marked as possible translocation", () => {
  const r = structuralEvidence(
    [read(0, { sa: "chr2,5001,+,50S50M,60,0;" })],
    region,
    opts,
  );
  assert.equal(r.events.length, 1);
  assert.equal(r.events[0].kind, "translocation-split");
  assert.equal(r.events[0].partner.pos, 5000);
  assert.equal(r.events[0].eligibleFragments, 1);
});
test("opposite-strand complementary split alignment is inversion-like", () => {
  const r = structuralEvidence(
    [read(0, { sa: "chr1,1501,-,50M50S,60,0;" })],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "inversion-split");
});
test("overlapping alternative mappings do not become inversion or translocation calls", () => {
  const r = structuralEvidence(
    [read(0, { sa: "chr2,5001,+,50M50S,60,0;" })],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "ambiguous-split");
});
test("supplementary records are checked without requiring SA and are not counted as separate fragments", () => {
  const a = read(0, { name: "split" }),
    b = read(1, {
      name: "split",
      start: 1500,
      end: 1550,
      cigar: "50S50M",
      flags: 2048,
    });
  const result = structuralEvidence([a, b], region, opts);
  assert.ok(result.events.some((e) => e.kind === "regional-split"));
  assert.ok(result.events.every((e) => e.fragments === 1));
});
test("mate-only interchromosomal links show unknown partner mapping quality", () => {
  const r = structuralEvidence(
    [read(0, { flags: 65, mate: 9000, mateChr: "chr8" })],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "translocation-pair");
  assert.equal(r.events[0].unknownPartnerQuality, 1);
  assert.equal(r.events[0].eligibleFragments, 0);
  assert.equal(r.events[0].localQualityFragments, 1);
});
test("same-strand mates are orientation clues, ordinary nearby FR mates are not flagged", () => {
  const normal = read(0, {
      flags: 97,
      mate: 1400,
      mateChr: "chr1",
      mateMq: 60,
    }),
    same = read(1, { flags: 65, mate: 1400, mateChr: "chr1", mateMq: 60 });
  const events = structuralEvidence([normal, same], region, opts).events;
  assert.equal(events.length, 1);
  assert.equal(events[0].kind, "inversion-pair");
  assert.match(events[0].notes[0], /libraries/);
});
test("distant intrachromosomal mates are shown without a translocation claim", () => {
  const r = structuralEvidence(
    [read(0, { flags: 97, mate: 90000, mateChr: "chr1", mateMq: 60 })],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "distant-pair");
});
test("support is deduplicated per fragment; duplicates and low MAPQ remain explicit weak evidence", () => {
  const a = read(0, { name: "a", sa: "chr2,5001,+,50S50M,60,0;" }),
    b = { ...a, id: 1, flags: 2048 },
    c = read(2, { flags: 1024, sa: a.sa }),
    d = read(3, { mq: 0, sa: a.sa });
  const e = structuralEvidence([a, b, c, d], region, opts).events.find(
    (e) => e.kind === "translocation-split",
  );
  assert.equal(e.fragments, 3);
  assert.equal(e.eligibleFragments, 1);
  assert.equal(e.duplicateFragments, 1);
});
test("fetched records outside the selected locus cannot seed unrelated signals", () => {
  const r = structuralEvidence(
    [read(0, { start: 5000, end: 5050, sa: "chr2,5001,+,50S50M,60,0;" })],
    region,
    opts,
  );
  assert.equal(r.events.length, 0);
});
test("structurally different alignments do not collapse into the same display family", () => {
  const reads = [
    read(0, { cigar: "100M", end: 1200, sa: "chr2,5001,+,50S50M,60,0;" }),
    read(1, { cigar: "100M", end: 1200 }),
  ];
  const r = phase(reads, "A".repeat(1000), 1000, opts);
  r.structural = structuralEvidence(r.reads, region, opts);
  assert.equal(buildFamilies(r, "A".repeat(1000), 1000, opts).length, 2);
});

test("hard-clipped supplementary coordinates use the original read orientation", () => {
  const r = structuralEvidence(
    [
      read(0, {
        cigar: "50M50H",
        seq: "A".repeat(50),
        sa: "chr1,1501,-,50M50H,60,0;",
      }),
    ],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "inversion-split");
});
test("large fractional overlap on short anchors stays ambiguous", () => {
  const r = structuralEvidence(
    [read(0, { cigar: "12M88S", end: 1112, sa: "chr2,5001,+,2S12M86S,60,0;" })],
    region,
    opts,
  );
  assert.equal(r.events[0].kind, "ambiguous-split");
});
test("mate MAPQ can be obtained from a record already available locally", () => {
  const a = read(0, { name: "pair", flags: 65, mate: 1300, mateChr: "chr1" }),
    b = read(1, {
      name: "pair",
      flags: 129,
      start: 1300,
      end: 1350,
      mate: 1100,
      mateChr: "chr1",
    });
  const events = structuralEvidence([a, b], region, opts).events;
  assert.ok(events.length);
  assert.ok(
    events.every(
      (e) => e.eligibleFragments === 1 && e.unknownPartnerQuality === 0,
    ),
  );
});
