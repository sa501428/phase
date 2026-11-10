import { test } from "node:test";
import assert from "node:assert/strict";
import { phase } from "../src/analysis.js";
import { buildFamilies, layoutFamilies } from "../src/families.js";
import { synthetic } from "../src/demo.js";
const opts = { mapq: 20, baseq: 20, support: 3, duplicates: false, sample: "" };
function model(reads, reference, origin = 0) {
  const result = phase(reads, reference, origin, opts);
  return { result, families: buildFamilies(result, reference, origin, opts) };
}
function read(id, name, start, seq, flags = 0) {
  return {
    id,
    name,
    start,
    end: start + seq.length,
    seq,
    qual: Array(seq.length).fill(35),
    cigar: seq.length + "M",
    flags,
    mq: 60,
    rg: "r",
    sample: "s",
    mate: -1,
  };
}
test("matching patterns collapse despite different endpoints without losing counts", () => {
  const d = synthetic(),
    { families, result } = model(d.reads, d.reference, 1000);
  assert.equal(families.length, 4);
  assert.deepEqual(
    families.map((f) => f.count).sort((a, b) => a - b),
    [8, 12, 12, 12],
  );
  assert.equal(
    families.reduce((n, f) => n + f.recordCount, 0),
    44,
  );
  assert.equal(result.fragmentCount, 44);
});
test("paired ends occupy one fragment row; expansion retains both mates of every fragment", () => {
  const reads = [
      read(0, "a", 0, "AAAA", 65),
      read(1, "a", 8, "AAAA", 129),
      read(2, "b", 1, "AAAA", 65),
      read(3, "b", 9, "AAAA", 129),
    ],
    { families } = model(reads, "A".repeat(20));
  assert.equal(families.length, 1);
  assert.equal(families[0].count, 2);
  assert.equal(families[0].recordCount, 4);
  let rows = layoutFamilies(families, reads, { start: 0, end: 20 });
  assert.equal(rows.length, 1);
  rows = layoutFamilies(families, reads, {
    start: 0,
    end: 20,
    expanded: new Set([0]),
  });
  assert.equal(rows.length, 3);
  assert.ok(
    rows.slice(1).every((r) => r.member.matePresent && r.reads.length === 2),
  );
});
test("duplicate-marked records collapse visually without changing analysis exclusion", () => {
  const reads = [read(0, "a", 0, "AAAA"), read(1, "b", 0, "AAAA", 1024)],
    { families, result } = model(reads, "AAAA");
  assert.equal(result.fragmentCount, 1);
  assert.equal(families.length, 1);
  assert.equal(families[0].duplicates, 1);
  assert.equal(families[0].count, 2);
  const rows = layoutFamilies(families, reads, {
    start: 0,
    end: 4,
    layer: "eligible",
  });
  assert.equal(rows[0].type, "fragment");
  assert.deepEqual(rows[0].member.readIds, [0]);
});
test("unknown candidate coverage and observed reference are not merged", () => {
  const reads = [
      read(0, "alt1", 0, "ATAAAA"),
      read(1, "alt2", 0, "ATAAAA"),
      read(2, "alt3", 0, "ATAAAA"),
      read(3, "ref", 0, "AAAAAA"),
      read(4, "unknown", 3, "AAA"),
    ],
    { families } = model(reads, "AAAAAA");
  assert.ok(
    !families.some(
      (f) =>
        f.members.some((m) => m.id === 3) && f.members.some((m) => m.id === 4),
    ),
  );
});
test("different rare SNPs, indels, clips and conflicts prevent a confident collapse", () => {
  const reads = [
      read(0, "ref", 0, "AAAAAA"),
      read(1, "rare", 0, "AACAAA"),
      { ...read(2, "ins", 0, "AAAGAAA"), cigar: "3M1I3M", end: 6 },
      { ...read(3, "clip", 0, "TTAAAAAA"), cigar: "2S6M", end: 6 },
    ],
    { families } = model(reads, "AAAAAA");
  assert.equal(families.length, 4);
});
test("expanded layout contains all fragments beyond the old 150-row limit", () => {
  const reads = Array.from({ length: 400 }, (_, i) =>
      read(i, "r" + i, 0, "AAAA"),
    ),
    { families } = model(reads, "AAAA");
  const rows = layoutFamilies(families, reads, {
    start: 0,
    end: 4,
    collapsed: false,
  });
  assert.equal(rows.length, 400);
  assert.equal(new Set(rows.flatMap((r) => r.member.readIds)).size, 400);
  assert.equal(rows.at(-1).member.id, 399);
});

test("disconnected reference-only coverage does not create a merged family", () => {
  const reads = [read(0, "left", 0, "AAAA"), read(1, "right", 8, "AAAA")];
  assert.equal(model(reads, "A".repeat(16)).families.length, 2);
});
test("identical low-quality duplicate records can collapse by complete raw identity", () => {
  const a = { ...read(0, "a", 0, "AACAAA"), qual: Array(6).fill(5) },
    b = { ...a, id: 1, name: "duplicate", flags: 1024 };
  const { families, result } = model([a, b], "AAAAAA");
  assert.equal(families.length, 1);
  assert.equal(families[0].count, 2);
  assert.equal(result.fragmentCount, 1);
});
test("an indel conflict between mates is not merged with consistent pairs", () => {
  const insert = (id, name, flags) => ({
    ...read(id, name, 0, "AAAGAAA", flags),
    cigar: "3M1I3M",
    end: 6,
  });
  const reads = [
    insert(0, "conflict", 65),
    read(1, "conflict", 0, "AAAAAA", 129),
    insert(2, "consistent", 65),
    insert(3, "consistent", 129),
  ];
  assert.equal(model(reads, "AAAAAA").families.length, 2);
});
