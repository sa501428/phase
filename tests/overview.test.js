import { test } from "node:test";
import assert from "node:assert/strict";
import { phase } from "../src/analysis.js";
import { buildFamilies, layoutFamilies } from "../src/families.js";
import { buildOverview } from "../src/overview.js";
import { consensusGeometry } from "../src/consensus.js";
import { synthetic } from "../src/demo.js";
const opts = {
  minVaf: 0,
  mapq: 20,
  baseq: 20,
  support: 3,
  duplicates: false,
  sample: "",
};
function read(id, start, seq, extra = {}) {
  return {
    id,
    name: "r" + id,
    start,
    end: start + seq.length,
    seq,
    qual: Array(seq.length).fill(35),
    cigar: seq.length + "M",
    flags: 0,
    mq: 60,
    rg: "r",
    sample: "s",
    mate: -1,
    ...extra,
  };
}
function model(reads, ref, origin = 0) {
  const r = phase(reads, ref, origin, opts);
  r.families = buildFamilies(r, ref, origin, opts);
  r.overview = buildOverview(
    r,
    { start: origin, end: origin + ref.length },
    opts,
  );
  return r;
}
test("rare SNP clashes split coverage rows while preserving every fragment", () => {
  const ref = "A".repeat(3200),
    reads = Array.from({ length: 3000 }, (_, i) =>
      read(i, i, "A".repeat(20) + "C" + "A".repeat(79)),
    ),
    r = model(reads, ref);
  assert.equal(r.families.length, 3000);
  assert.ok(r.overview.length > 1);
  assert.equal(
    r.overview.reduce((n, g) => n + g.count, 0),
    3000,
  );
  const rows = layoutFamilies(r.overview, r.reads, {
    start: 0,
    end: ref.length,
  });
  assert.equal(rows.length, r.overview.length);
  const expanded = layoutFamilies(r.overview, r.reads, {
    start: 0,
    end: ref.length,
    expanded: new Set(r.overview.map((g) => g.id)),
  });
  assert.equal(expanded.length, 3000 + r.overview.length);
  assert.equal(
    new Set(
      expanded
        .filter((row) => row.type === "fragment")
        .flatMap((row) => row.member.readIds),
    ).size,
    3000,
  );
});
test("partial compatible reads join a uniquely supported haplotype despite errors and clipping", () => {
  const ref = "A".repeat(40),
    reads = [];
  for (let i = 0; i < 4; i++) {
    let seq = ref.split("");
    seq[10] = "C";
    seq[20] = "G";
    reads.push(read(reads.length, 0, seq.join("")));
  }
  for (let i = 0; i < 4; i++) reads.push(read(reads.length, 0, ref));
  // One C allele uniquely identifies the C/G pattern; position 20 is uncovered.
  reads.push(read(8, 5, "AAAAACAAAA"));
  reads.push(read(9, 5, "TTAAAAACAAAA", { cigar: "2S10M", end: 15 }));
  const r = model(reads, ref),
    group = r.overview.find(
      (g) => g.category === "haplotype" && g.pattern[0] === "C",
    );
  assert.ok(group);
  assert.equal(group.count, 6);
  assert.equal(r.overview.filter((g) => g.category === "haplotype").length, 2);
  assert.ok(
    group.members.some((m) => m.id === 8) &&
      group.members.some((m) => m.id === 9),
  );
});
test("ambiguous partials stay in a single bucket instead of inventing haplotypes", () => {
  const ref = "A".repeat(40),
    reads = [];
  for (const base of ["G", "T"])
    for (let i = 0; i < 3; i++) {
      const seq = ref.split("");
      seq[10] = "C";
      seq[20] = base;
      reads.push(read(reads.length, 0, seq.join("")));
    }
  for (let i = 0; i < 100; i++) reads.push(read(reads.length, 5, "AAAAACAAAA"));
  const r = model(reads, ref);
  assert.equal(r.overview.filter((g) => g.category === "haplotype").length, 2);
  assert.equal(r.overview.find((g) => g.category === "ambiguous").count, 100);
  assert.equal(r.overview.length, 3);
});
test("duplicates are absorbed without contributing to analysis support", () => {
  const ref = "A".repeat(20),
    reads = [];
  for (let i = 0; i < 3; i++)
    reads.push(read(i, 0, "A".repeat(10) + "C" + "A".repeat(9)));
  reads.push({ ...reads[0], id: 3, name: "dup", flags: 1024 });
  const r = model(reads, ref);
  assert.equal(r.fragmentCount, 3);
  assert.equal(r.overview.length, 1);
  assert.equal(r.overview[0].count, 4);
  assert.equal(r.overview[0].duplicates, 1);
});
test("summary consensus flags contradictions without a majority override and counts overlapping mates once", () => {
  const ref = "A".repeat(20),
    reads = Array.from({ length: 10 }, (_, i) => read(i, 0, "A".repeat(20))),
    byID = new Map(reads.map((r) => [r.id, r]));
  reads[0].seq = "CA" + "A".repeat(18);
  for (let i = 0; i < 4; i++)
    reads[i].seq = reads[i].seq.slice(0, 10) + "G" + reads[i].seq.slice(11);
  const glyph = consensusGeometry(
    reads.map((r) => ({ readIds: [r.id] })),
    byID,
    ref,
    0,
    opts,
  );
  assert.equal(glyph.mismatches.length, 0); // 4/10 is not a confident consensus
  for (let i = 4; i < 8; i++)
    reads[i].seq = reads[i].seq.slice(0, 10) + "G" + reads[i].seq.slice(11);
  const supported = consensusGeometry(
    reads.map((r) => ({ readIds: [r.id] })),
    byID,
    ref,
    0,
    opts,
  );
  assert.equal(supported.mismatches.length, 0);
  assert.deepEqual(
    supported.conflicts.map((c) => c.pos).sort((a, b) => a - b),
    [0, 10],
  );
  const one = consensusGeometry([{ readIds: [0, 1] }], byID, ref, 0, opts);
  assert.equal(one.mismatches.find((m) => m.pos === 10).support, 1);
  assert.ok(one.conflicts.some((c) => c.pos === 0));
});
test("disconnected phase blocks stay separate and every input record remains accessible", () => {
  const d = synthetic(),
    r = model(d.reads, d.reference, 1000);
  assert.equal(r.overview.length, 4);
  assert.deepEqual([...new Set(r.overview.map((g) => g.blockID))], [1, 2]);
  assert.equal(
    new Set(r.overview.flatMap((g) => g.members.flatMap((m) => m.readIds)))
      .size,
    44,
  );
});

test("one thousand clipping families are absorbed into one supported haplotype group", () => {
  const ref = "A".repeat(80),
    reads = Array.from({ length: 1000 }, (_, i) => {
      const bases = ref.split("");
      bases[40] = "C";
      return read(i, 0, "T".repeat(i) + bases.join(""), {
        cigar: (i ? i + "S" : "") + "80M",
        end: 80,
      });
    }),
    r = model(reads, ref);
  assert.equal(r.families.length, 1000);
  assert.equal(r.overview.length, 1);
  assert.equal(r.overview[0].category, "haplotype");
  assert.equal(r.overview[0].count, 1000);
});
test("supported distinct structural signals stay separate even without SNP candidates", () => {
  const ref = "A".repeat(40),
    reads = Array.from({ length: 9 }, (_, i) => read(i, 0, ref)),
    r = phase(reads, ref, 0, opts);
  for (let i = 0; i < 6; i++) r.reads[i].structuralLinks = [i < 3 ? 1 : 2];
  r.structural = {
    events: [
      { id: 1, kind: "translocation-split", eligibleFragments: 3 },
      { id: 2, kind: "inversion-split", eligibleFragments: 3 },
    ],
  };
  r.families = buildFamilies(r, ref, 0, opts);
  const overview = buildOverview(r, { start: 0, end: 40 }, opts);
  assert.equal(overview.length, 3);
  assert.equal(
    new Set(overview.flatMap((g) => g.members.flatMap((m) => m.readIds))).size,
    9,
  );
});
