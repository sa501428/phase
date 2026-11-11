import { test } from "node:test";
import assert from "node:assert/strict";
import { applyVafFilter, hiddenBase, hiddenEvent } from "../src/vaf-filter.js";
import { ManualWorkspace } from "../v3/model.js";
import { phase, normalize } from "../src/analysis.js";
import { buildFamilies } from "../src/families.js";
import { rawGeometry } from "../src/alignment-geometry.js";
import { assemble } from "../v2/assembly.js";
import { somaticDemo } from "../v2/somatic-demo.js";
import { filterRequestedVariants } from "../src/variant-input.js";
import { structuralEvidence } from "../src/structural.js";
const reference = "ACGT".repeat(25),
  region = { chr: "chr1", start: 1000, end: 1100 };
const opts = { mapq: 20, baseq: 20, support: 1, duplicates: false };
function read(id, extra = {}) {
  return {
    id,
    name: "r" + id,
    chr: "chr1",
    sample: "s",
    rg: "r",
    start: 1000,
    end: 1100,
    seq: reference,
    cigar: "100M",
    qual: Array(100).fill(35),
    mq: 60,
    flags: 0,
    mate: -1,
    ...extra,
  };
}
function mixture(n = 51) {
  const reads = Array.from({ length: n }, (_, i) => read(i));
  reads.at(-1).seq = reference.slice(0, 30) + "A" + reference.slice(31);
  return reads;
}
function workspace(reads) {
  return new ManualWorkspace({ reads, region, reference, dataset: "vaf-test" });
}
test("default cutoff hides <2%, retains exactly 2%, and zero restores observations without changing sequence", () => {
  const at = mixture(50),
    below = mixture(51);
  applyVafFilter(at, region, reference, opts);
  assert.equal(hiddenBase(at.at(-1), 30), false);
  const original = below.at(-1).seq;
  assert.equal(applyVafFilter(below, region, reference, opts).suppressed, 1);
  assert.equal(hiddenBase(below.at(-1), 30), true);
  assert.equal(below.length, 51);
  assert.equal(below.at(-1).seq, original);
  applyVafFilter(below, region, reference, { ...opts, minVaf: 0 });
  assert.equal(hiddenBase(below.at(-1), 30), false);
  applyVafFilter(at, region, reference, { ...opts, minVaf: 0.03 });
  assert.equal(hiddenBase(at.at(-1), 30), true);
});
test("physical pairs count once; off-contig, low quality and duplicate observations do not dilute local VAF", () => {
  const reads = mixture(50);
  reads.push(
    ...reads.slice(0, -1).map((r) => ({ ...r, id: r.id + 50, flags: 128 })),
  );
  reads.push(
    ...Array.from({ length: 80 }, (_, i) => read(200 + i, { chr: "chr2" })),
  );
  reads.push(...Array.from({ length: 80 }, (_, i) => read(300 + i, { mq: 0 })));
  reads.push(
    ...Array.from({ length: 80 }, (_, i) => read(400 + i, { flags: 1024 })),
  );
  reads.push(
    ...Array.from({ length: 80 }, (_, i) =>
      read(500 + i, { qual: Array(100).fill(0) }),
    ),
  );
  applyVafFilter(reads, region, reference, opts);
  assert.equal(
    hiddenBase(reads[49], 30),
    false,
    "1/50 pair VAF stays at boundary",
  );
});
test("per-site coverage, not all fragments in the locus, determines suppression", () => {
  const reads = mixture(10);
  reads.push(
    ...Array.from({ length: 100 }, (_, i) =>
      read(20 + i, {
        start: 1060,
        end: 1100,
        seq: reference.slice(60),
        cigar: "40M",
        qual: Array(40).fill(35),
      }),
    ),
  );
  applyVafFilter(reads, region, reference, opts);
  assert.equal(hiddenBase(reads[9], 30), false);
});
test("rare alternate is unknown for manual matching, never reference; restored allele vetoes merge", () => {
  const m = workspace(mixture());
  assert.equal(m.settings.minVaf, 0.02);
  assert.equal(m.evidence(50).bases.has(1030), false);
  assert.equal(
    rawGeometry([m.reads[50]], reference, 1000).mismatches.length,
    0,
  );
  assert.equal(m.evidence(0).bases.get(1030).has("G"), true);
  m.selection = new Set([0, 50]);
  m.merge("rare is ignored");
  assert.equal(m.groups.length, 1);
  m.setSettings({ ...m.settings, minVaf: 0 });
  assert.equal(
    m.groupSummary(m.groups[0]).bases.find((b) => b.pos === 1030).conflict,
    true,
  );
  const restored = workspace(mixture());
  restored.setSettings({ ...restored.settings, minVaf: 0 });
  restored.selection = new Set([0, 50]);
  assert.throws(() => restored.merge(), /contradictory/);
});
test("common multiallelic SNPs still clash when another alternate at that position is suppressed", () => {
  const reads = mixture(100);
  for (let i = 0; i < 10; i++)
    reads[i].seq = reference.slice(0, 30) + "T" + reference.slice(31);
  const m = workspace(reads);
  assert.equal(m.evidence(99).bases.has(1030), false);
  assert.equal(m.evidence(0).bases.get(1030).has("T"), true);
  m.selection = new Set([0, 10]);
  assert.throws(() => m.merge(), /contradictory/);
});
for (const op of ["I", "D"])
  test(`rare ${op} is hidden from manual evidence and rendering and cannot clash`, () => {
    const reads = mixture();
    reads.at(-1).seq = reference;
    const r = reads.at(-1);
    if (op === "I") {
      r.seq = reference.slice(0, 30) + "TT" + reference.slice(30);
      r.cigar = "30M2I70M";
    } else {
      r.seq = reference.slice(0, 30) + reference.slice(32);
      r.cigar = "30M2D68M";
    }
    r.qual = Array(r.seq.length).fill(35);
    const m = workspace(reads);
    assert.equal(hiddenEvent(r, op, 1030), true);
    assert.equal(rawGeometry([r], reference, 1000).events.length, 0);
    assert.equal(
      m.evidence(50).insertions.get(1030)?.has("TT") || false,
      false,
    );
    assert.equal(m.evidence(50).bases.get(1030)?.has("-") || false, false);
    m.selection = new Set([0, 50]);
    m.merge();
    m.setSettings({ ...m.settings, minVaf: 0 });
    assert.equal(hiddenEvent(r, op, 1030), false);
    assert.equal(
      m
        .groupSummary(m.groups[0])
        [op === "I" ? "insertions" : "bases"].some((c) => c.conflict),
      true,
    );
  });
test("one-sided insertion coverage is unknown and cannot dilute an observed allele", () => {
  const reads = [
    read(0, {
      seq: reference.slice(0, 30) + "TT" + reference.slice(30),
      cigar: "30M2I70M",
      qual: Array(102).fill(35),
    }),
  ];
  reads.push(
    ...Array.from({ length: 80 }, (_, i) =>
      read(i + 1, {
        start: 1030,
        seq: reference.slice(30),
        cigar: "70M",
        qual: Array(70).fill(35),
      }),
    ),
  );
  applyVafFilter(reads, region, reference, opts);
  assert.equal(hiddenEvent(reads[0], "I", 1030), false);
});
test("normalized insertion alleles in repeats do not contradict their own reference steps", () => {
  const ref = "A".repeat(100),
    reads = Array.from({ length: 51 }, (_, i) => read(i, { seq: ref }));
  reads.at(-1).seq = "A".repeat(101);
  reads.at(-1).cigar = "30M1I70M";
  reads.at(-1).qual = Array(101).fill(35);
  // Absence can be established at the equivalent raw anchor even when the
  // normalized event is at the target edge. Its own steps are never REF votes.
  applyVafFilter(reads, region, ref, opts);
  assert.equal(hiddenEvent(reads.at(-1), "I", 1030), true);
});
test("v1 suppresses a supplied rare variant too, and does not retain a family split or SNP glyph", () => {
  const reads = mixture(),
    supplied = [normalize(1030, "G", "A", reference, 1000)];
  const r = phase(reads, reference, 1000, opts, supplied);
  assert.equal(r.variants.length, 0);
  const requested = filterRequestedVariants(
    [{ variant: supplied[0], input: "chr1:1031G>A" }],
    r,
  );
  assert.equal(requested[0].variant, null);
  assert.match(requested[0].error, /Suppressed below minimum VAF/);
  assert.equal(buildFamilies(r, reference, 1000, opts).length, 1);
  const restored = phase(
    reads,
    reference,
    1000,
    { ...opts, minVaf: 0 },
    supplied,
  );
  assert.equal(restored.variants.length, 1);
  assert.equal(buildFamilies(restored, reference, 1000, opts).length, 2);
});
test("assembly default suppresses a <2% somatic private branch; changing cutoff restores ten lineages", () => {
  const d = somaticDemo(),
    r = assemble(d.reads, d.region, d.reference);
  assert.equal(r.settings.minVaf, 0.02);
  assert.equal(r.blocks[0].groups.filter((g) => g.supported).length, 9);
  assert.equal(r.stats.records, d.reads.length);
  assert.ok(
    r.fragments
      .filter((f) => f.name.startsWith("lineage-10-"))
      .every((f) => !new Map(f.calls).has(1450)),
  );
  const all = assemble(d.reads, d.region, d.reference, { minVaf: 0 });
  assert.equal(all.blocks[0].groups.filter((g) => g.supported).length, 10);
});
test("rare structural clues are omitted from annotations and family fingerprints and can be restored", () => {
  const reads = Array.from({ length: 51 }, (_, i) =>
    read(i, { flags: 65, mate: 9000, mateChr: "chr1", mateMq: 60 }),
  );
  reads.at(-1).mateChr = "chr2";
  const result = structuralEvidence(reads, region, opts);
  assert.ok(!result.events.some((e) => e.kind === "translocation-pair"));
  assert.deepEqual(reads.at(-1).structuralLinks, []);
  assert.ok(
    structuralEvidence(reads, region, { ...opts, minVaf: 0 }).events.some(
      (e) => e.kind === "translocation-pair",
    ),
  );
});
test("invalid VAF values are rejected", () => {
  for (const v of [-1, 1.1, NaN, Infinity])
    assert.throws(
      () => applyVafFilter([], region, reference, { minVaf: v }),
      /VAF/,
    );
});

test("zero-quality rare observations are suppressed when callable fragments establish zero eligible ALT support", () => {
  const reads = mixture(3);
  reads.at(-1).qual = Array(100).fill(0);
  const r = phase(reads, reference, 1000, opts);
  assert.equal(hiddenBase(reads.at(-1), 30), true);
  assert.equal(buildFamilies(r, reference, 1000, opts).length, 1);
});

test("assembly does not render rare soft-clip-derived allele branches, and lowering cutoff restores them", async () => {
  const { demo } = await import("../v2/demo.js");
  const { reverseComplement } = await import("../v2/assembly.js");
  const ref = demo().reference.slice(0, 100),
    alt = ref.slice(0, 20) + (ref[20] === "A" ? "T" : "A") + ref.slice(21),
    marker = alt.slice(13, 28);
  const reads = Array.from({ length: 100 }, (_, id) =>
    read(id, {
      chr: "synthetic",
      start: 1050,
      end: 1100,
      cigar: "50S50M",
      seq: id === 0 ? alt : ref,
    }),
  );
  const target = { chr: "synthetic", start: 1000, end: 1100 };
  const hasBranch = (r) =>
    r.unitigs.some(
      (u) =>
        u.seq.includes(marker) || u.seq.includes(reverseComplement(marker)),
    );
  const filtered = assemble(reads, target, ref, { k: 15 });
  assert.equal(hasBranch(filtered), false);
  assert.equal(filtered.bubbles.length, 0);
  assert.equal(
    hasBranch(assemble(reads, target, ref, { k: 15, minVaf: 0 })),
    true,
  );
});
test("assembly retains inspectable reads when suppression leaves no k-mer path", () => {
  const ref = "A".repeat(100),
    reads = Array.from({ length: 4 }, (_, id) =>
      read(id, {
        seq: Array.from({ length: 100 }, (_, q) => "ACGT"[(q + id) % 4]).join(
          "",
        ),
      }),
    );
  const filtered = assemble(reads, region, ref, { k: 15, minVaf: 0.3 });
  assert.equal(filtered.unitigs.length, 0);
  assert.equal(filtered.reads.length, 4);
  assert.ok(
    assemble(reads, region, ref, { k: 15, minVaf: 0 }).unitigs.length > 0,
  );
});
