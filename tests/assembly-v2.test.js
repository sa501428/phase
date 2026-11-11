import { test } from "node:test";
import assert from "node:assert/strict";
import { assemble, exportGFA, reverseComplement } from "../v2/assembly.js";
import { demo } from "../v2/demo.js";
import { queryVariants } from "../v2/queries.js";
const build = (reference = true, opts = {}) => {
  const d = demo();
  return assemble(d.reads, d.region, reference ? d.reference : null, opts);
};
test("assembly resolves three linked sequence patterns and two deletion alleles without a reference", () => {
  const r = build(false);
  assert.equal(r.stats.records, 104);
  assert.equal(r.stats.fragments, 52);
  assert.equal(r.bubbles.length, 4);
  assert.deepEqual(
    r.blocks.map((b) => b.groups.map((g) => g.fragments.length)),
    [
      [12, 12, 12],
      [8, 8],
    ],
  );
  assert.ok(r.bubbles.every((b) => !b.anchored));
  assert.equal(r.blocks[0].total, 36);
  assert.equal(r.blocks[1].total, 16);
});
test("optional reference annotates existing graph and never supplies extra sequence", () => {
  const a = build(false),
    b = build(true);
  assert.deepEqual(
    a.unitigs.map((u) => u.seq),
    b.unitigs.map((u) => u.seq),
  );
  assert.deepEqual(a.links, b.links);
  assert.deepEqual(a.blocks, b.blocks);
  assert.ok(b.bubbles.every((b) => b.anchored));
  assert.ok(b.bubbles.every((b) => b.alleles.some((a) => a.isReference)));
});
test("unitig export preserves oriented overlaps without doubling reverse-complement segments", () => {
  const r = build();
  for (const e of r.links) {
    const a =
        e.fromOrientation === "+"
          ? r.unitigs[e.from].seq
          : reverseComplement(r.unitigs[e.from].seq),
      b =
        e.toOrientation === "+"
          ? r.unitigs[e.to].seq
          : reverseComplement(r.unitigs[e.to].seq);
    assert.equal(a.slice(-(r.settings.k - 1)), b.slice(0, r.settings.k - 1));
  }
  const canonical = r.unitigs.map(
    (u) => [u.seq, reverseComplement(u.seq)].sort()[0],
  );
  assert.equal(new Set(canonical).size, canonical.length);
  assert.match(exportGFA(r), /^H\tVN:Z:1.0\nS\tu1\t[ACGT]+/);
  assert.equal(
    exportGFA(r)
      .split("\n")
      .filter((s) => s.startsWith("L\t")).length,
    r.links.length,
  );
});
test("one private contradictory SNP cannot enter a group containing the opposite base", () => {
  const d = demo(),
    r = d.reads[0],
    q = 1250 - r.start;
  r.seq =
    r.seq.slice(0, q) + (r.seq[q] === "A" ? "T" : "A") + r.seq.slice(q + 1);
  const result = assemble(d.reads, d.region, null, { minVaf: 0 });
  for (const b of result.blocks)
    for (const g of b.groups) {
      const calls = new Map();
      for (const id of g.fragments)
        for (const [pos, base] of result.fragments[id].calls) {
          assert.ok(
            !calls.has(pos) || calls.get(pos) === base,
            `Contradictory ${pos}`,
          );
          calls.set(pos, base);
        }
    }
  assert.ok(
    result.fragments.find((f) => f.name === r.name).conflict,
    "the mate still has the contradictory base",
  );
});
test("singleton conflict veto holds when the mate does not cover the site", () => {
  const d = demo(),
    r = d.reads[0],
    q = 45;
  r.seq =
    r.seq.slice(0, q) + (r.seq[q] === "A" ? "T" : "A") + r.seq.slice(q + 1);
  const result = assemble(d.reads, d.region, null);
  for (const b of result.blocks)
    for (const g of b.groups) {
      const seen = new Map();
      for (const id of g.fragments)
        for (const [pos, base] of result.fragments[id].calls) {
          assert.ok(!seen.has(pos) || seen.get(pos) === base);
          seen.set(pos, base);
        }
    }
  assert.ok(
    result.blocks.some((b) =>
      b.groups.some((g) => g.fragments.length === 1 && !g.supported),
    ),
  );
});
test("missing sequence, missing CIGAR, unmapped and inconsistent records are quarantined without aborting", () => {
  const d = demo(),
    base = d.reads[0];
  for (const [i, change] of [
    { seq: "" },
    { cigar: "" },
    { flags: 4, cigar: "" },
    { seq: base.seq.slice(1) },
  ].entries())
    d.reads.push({ ...base, ...change, id: 104 + i, name: "bad-" + i });
  const r = assemble(d.reads, d.region, null);
  assert.equal(r.stats.records, 108);
  assert.equal(r.stats.eligibleRecords, 104);
  assert.equal(r.stats.issues["missing-sequence"], 1);
  assert.equal(r.stats.issues["inconsistent-sequence-length"], 1);
  assert.equal(r.blocks[0].groups.length, 3);
});
test("duplicate flags change voting only when explicitly included; input records stay retained", () => {
  const d = demo();
  d.reads[0].flags |= 1024;
  d.reads[1].flags |= 1024;
  const a = assemble(structuredClone(d.reads), d.region, null),
    b = assemble(structuredClone(d.reads), d.region, null, {
      duplicates: true,
    });
  assert.equal(a.stats.records, b.stats.records);
  assert.equal(a.stats.fragments, 51);
  assert.equal(b.stats.fragments, 52);
});
test("weak phase bridge does not manufacture a long haplotype block", () => {
  const d = demo(),
    reference = d.reference,
    region = { chr: "synthetic", start: 1000, end: 1800 },
    reads = [];
  const add = (name, start, alt, flags = 0) => {
    let seq = reference.slice(start - 1000, start - 1000 + 100);
    const pos = start < 1200 ? 1120 : 1500;
    const q = pos - start;
    if (alt && q >= 0 && q < 100)
      seq = seq.slice(0, q) + (seq[q] === "A" ? "T" : "A") + seq.slice(q + 1);
    reads.push({
      id: reads.length,
      name,
      start,
      end: start + 100,
      seq,
      cigar: "100M",
      qual: Array(100).fill(35),
      mq: 60,
      flags,
    });
  };
  for (const start of [1070, 1450])
    for (let a = 0; a < 2; a++)
      for (let n = 0; n < 5; n++) add(`${start}-${a}-${n}`, start, a);
  add("bridge", 1070, true, 65);
  add("bridge", 1450, true, 129);
  const r = assemble(reads, region, null);
  assert.equal(r.bubbles.length, 2);
  assert.equal(r.blocks.length, 2);
  assert.equal(r.phaseLinks[0].support, 1);
  assert.equal(r.phaseLinks[0].linked, false);
});
test("reference-free queries count paired fragments once and never assume diploid trans", () => {
  const r = build(false),
    d = demo(),
    ref1 = d.reference[120],
    ref2 = d.reference[160],
    alt1 = ref1 === "A" ? "T" : "A",
    alt2 = ref2 === "C" ? "G" : "C";
  const q = queryVariants(
    `synthetic:1121${ref1}>${alt1}\nsynthetic:1161${ref2}>${alt2}`,
    r,
  );
  assert.equal(q.pairs[0].both, 24);
  assert.equal(q.pairs[0].neither, 12);
  assert.match(q.pairs[0].status, /Cis/);
  assert.ok(q.entries.every((e) => e.unverified));
});
test("reference-backed queries retain HGVS resolver and same physical linkage counts", () => {
  const r = build(),
    d = demo(),
    ref1 = d.reference[120],
    ref2 = d.reference[160],
    alt1 = ref1 === "A" ? "T" : "A",
    alt2 = ref2 === "C" ? "G" : "C";
  const q = queryVariants(
    `synthetic:g.1121${ref1}>${alt1}\nsynthetic:g.1161${ref2}>${alt2}`,
    r,
  );
  assert.equal(q.pairs[0].both, 24);
  assert.equal(q.pairs[0].neither, 12);
});
test("a single invariant graph path is not presented as a measured clone count", () => {
  const d = demo();
  const reads = d.reads.filter((r) => r.name.startsWith("haplotype-1"));
  const r = assemble(reads, d.region, null);
  assert.equal(r.bubbles.length, 0);
  assert.equal(r.blocks.length, 0);
  assert.equal(r.stats.unthreaded, 12);
});
test("repeat cycles remain in graph rather than being counted as extra clones", () => {
  const seq = "ACG".repeat(60),
    reads = Array.from({ length: 5 }, (_, id) => ({
      id,
      name: "repeat-" + id,
      start: 1000,
      end: 1180,
      seq,
      cigar: "180M",
      qual: Array(180).fill(35),
      flags: 0,
      mq: 60,
    }));
  const r = assemble(reads, { chr: "r", start: 1000, end: 1200 }, null, {
    k: 21,
  });
  assert.ok(r.unitigs.some((u) => u.cyclic));
  assert.ok(r.links.some((e) => e.from === e.to));
  assert.equal(r.blocks.length, 0);
});

test("ten somatic lineages of unequal abundance resolve without a clone-count or ploidy input", async () => {
  const { somaticDemo } = await import("../v2/somatic-demo.js"),
    d = somaticDemo(),
    r = assemble(d.reads, d.region, null, { minVaf: 0 });
  assert.equal(r.blocks.length, 1);
  assert.equal(r.blocks[0].groups.filter((g) => g.supported).length, 10);
  assert.deepEqual(
    r.blocks[0].groups.map((g) => g.fragments.length).sort((a, b) => a - b),
    d.truth.abundances.slice().sort((a, b) => a - b),
  );
  assert.equal(r.stats.fragments, 174);
  for (const g of r.blocks[0].groups) {
    const names = new Set(
      g.fragments.map((id) => r.fragments[id].name.match(/^lineage-\d+/)[0]),
    );
    assert.equal(names.size, 1, "never mix sibling or nested lineages");
  }
  const withRef = assemble(structuredClone(d.reads), d.region, d.reference, {
    minVaf: 0,
  });
  assert.deepEqual(r.blocks, withRef.blocks);
});
test("conserved founder-only observations stay shared among related somatic subclones", async () => {
  const { somaticDemo } = await import("../v2/somatic-demo.js"),
    d = somaticDemo({ partials: true }),
    r = assemble(d.reads, d.region, null, { minVaf: 0 }),
    b = r.blocks[0];
  assert.equal(b.groups.filter((g) => g.supported).length, 10);
  assert.equal(b.shared.length, 5);
  assert.equal(
    b.groups.reduce((n, g) => n + g.fragments.length, 0),
    174,
  );
  assert.equal(b.total, 179);
  assert.ok(b.shared.every((s) => s.groups.length === 5));
});
test("nested mutations retain cis linkage; sibling private mutations do not become diploid trans", async () => {
  const { somaticDemo } = await import("../v2/somatic-demo.js"),
    d = somaticDemo(),
    r = assemble(d.reads, d.region, null),
    entry = (p) =>
      `synthetic:${p + 1}${d.reference[p - 1000]}>${d.truth.alt(p)}`;
  const nested = queryVariants([1050, 1130, 1410].map(entry).join("\n"), r);
  assert.equal(
    nested.pairs.find((p) => p.left === entry(1050) && p.right === entry(1130))
      .both,
    30,
  );
  assert.equal(
    nested.pairs.find((p) => p.left === entry(1130) && p.right === entry(1410))
      .both,
    6,
  );
  const sibling = queryVariants([1130, 1170].map(entry).join("\n"), r).pairs[0];
  assert.equal(sibling.both, 0);
  assert.equal(sibling.leftOnly, 30);
  assert.equal(sibling.rightOnly, 18);
  assert.match(sibling.status, /trans is not established/);
});
test("off-target mate sequences do not assemble another region, while both records remain inspectable", () => {
  const d = demo(),
    baseline = assemble(structuredClone(d.reads), d.region, null);
  const r = d.reads[0];
  d.reads.push({
    ...r,
    id: d.reads.length,
    start: 5000,
    end: 5220,
    flags: 147,
    mate: r.start,
  });
  const assembled = assemble(d.reads, d.region, null);
  assert.deepEqual(
    assembled.unitigs.map((u) => u.seq),
    baseline.unitigs.map((u) => u.seq),
  );
  assert.equal(assembled.stats.issues["outside-target"], 1);
  assert.ok(assembled.fragments[0].reads.includes(104));
});
test("assembly uses only the requested interval, rather than the entire aligned read or BAM", () => {
  const d = demo(),
    region = { ...d.region, start: 1060, end: 1110 },
    r = assemble(d.reads, region, null);
  assert.equal(r.bubbles.length, 0);
  for (const u of r.unitigs)
    assert.ok(
      d.reference.slice(60, 110).includes(u.seq) ||
        d.reference.slice(60, 110).includes(reverseComplement(u.seq)),
    );
  assert.ok(r.stats.issues["outside-target"] > 0);
});

test("a rare high-quality error in a ten-lineage mixture stays weak instead of contaminating a supported group", async () => {
  const { somaticDemo } = await import("../v2/somatic-demo.js"),
    d = somaticDemo(),
    read = d.reads[0],
    q = 1100 - read.start;
  read.seq =
    read.seq.slice(0, q) +
    (read.seq[q] === "A" ? "T" : "A") +
    read.seq.slice(q + 1);
  const r = assemble(d.reads, d.region, null, { minVaf: 0 }),
    groups = r.blocks.flatMap((b) => b.groups);
  assert.equal(groups.filter((g) => g.supported).length, 10);
  assert.ok(groups.some((g) => !g.supported && g.fragments.length === 1));
  for (const g of groups) {
    const seen = new Map();
    for (const id of g.fragments)
      for (const [p, b] of r.fragments[id].calls) {
        assert.ok(!seen.has(p) || seen.get(p) === b);
        seen.set(p, b);
      }
  }
});
