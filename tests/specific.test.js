import test from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { parseVariants, resolveContigs } from "../specific/variants.js";
import {
  prepareRead,
  callVariant,
  makeCollector,
  summarize,
  pairMetrics,
  windowsForVariants,
} from "../specific/phasing.js";
import {
  discoverHelpers,
  helperGraph,
  findHelperPath,
  analyzeReads,
  helperWindows,
} from "../specific/helpers.js";
const opts = {
  mapq: 20,
  baseq: 20,
  duplicates: false,
  sample: "test",
  helpers: false,
  helperSupport: 5,
  helperFraction: 0.02,
  helperGap: 50000,
  linkAgreement: 0.95,
};
const read = (fields = {}) => ({
  name: "read",
  start: 100,
  seq: "A".repeat(100),
  qual: Array(100).fill(35),
  cigar: "100M",
  flags: 0,
  mq: 60,
  chr: "chr1",
  rg: "rg",
  sample: "test",
  ...fields,
});
const parse = (text) => parseVariants(text, { assembly: "hg38" });
const v = (pos, alt = "C", id = 0) => ({
  id,
  chr: "chr1",
  pos,
  ref: "A",
  alt,
  refLength: 1,
  label: `chr1:${pos + 1}A>${alt}`,
});
function paired(name, x, y, a, b, extras = {}) {
  return [
    read({
      name,
      start: x,
      seq: a,
      qual: Array(a.length).fill(35),
      cigar: `${a.length}M`,
      flags: 99,
      mate: y,
      mateChr: "chr1",
      ...extras,
    }),
    read({
      name,
      start: y,
      seq: b,
      qual: Array(b.length).fill(35),
      cigar: `${b.length}M`,
      flags: 147,
      mate: x,
      mateChr: "chr1",
      ...extras,
    }),
  ];
}
function sequence(start, length, changes) {
  const s = Array(length).fill("A");
  for (const [p, b] of changes)
    if (p >= start && p < start + length) s[p - start] = b;
  return s.join("");
}
export function helperReads() {
  const reads = [];
  for (const [block, x, y, changes] of [
    [
      0,
      100,
      180,
      [
        [120, "C"],
        [200, "G"],
      ],
    ],
    [
      1,
      180,
      580,
      [
        [200, "G"],
        [600, "T"],
      ],
    ],
    [
      2,
      580,
      680,
      [
        [600, "T"],
        [700, "C"],
      ],
    ],
  ]) {
    for (let i = 0; i < 12; i++) {
      const xa = x + i,
        ya = y + i,
        cs = i < 6 ? changes : [];
      reads.push(
        ...paired(
          `block${block}-${i}`,
          xa,
          ya,
          sequence(xa, 50, cs),
          sequence(ya, 50, cs),
        ),
      );
    }
  }
  return reads;
}

test("specific parser accepts genomic tables, VCF, HGVS and coalesces anchored alleles", () => {
  const result = parse(
    "chr1:101A>G\n1:101:A:G\nchr1 101 A G\nchr1\t101\t.\tA\tG,T\t.\tPASS\t.\nchr1:g.105_106delAA\nchr1:g.106_107insCG\nchr1:105:AAA:A",
  );
  assert.equal(
    result.entries.filter((e) => e.error).length,
    0,
    JSON.stringify(result.entries),
  );
  assert.equal(result.variants.length, 4); // G,T, deletion, insertion; the VCF deletion coalesces
  assert.equal(result.entries.at(-1).duplicate, true);
  assert.equal(result.entries[1].duplicate, true);
  assert.equal(result.variants[2].pos, 104);
  assert.equal(result.variants[3].pos, 106);
});
test("headed CSV/MAF columns and split ALT retain exact input line numbers", () => {
  const r = parse(
    '\n# comment\nHugo_Symbol,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2\nTP53,17,100,A,"C,G"\nTP53,17,101,A,T',
  );
  assert.equal(r.variants.length, 3);
  assert.ok(r.entries.every((e) => !e.error));
  assert.equal(r.entries[0].line, 4);
});
test("unlimited variant list replaces the old 32-variant restriction and preserves bad lines", () => {
  const r = parse(
    Array.from({ length: 200 }, (_, i) => `chr1 ${101 + i} A C`).join("\n") +
      "\np.R1G\nrs123",
  );
  assert.equal(r.variants.length, 200);
  assert.deepEqual(
    r.entries.filter((e) => e.error).map((e) => e.line),
    [201, 202],
  );
});
test("coding HGVS uses strand-aware transcript context and rejects ambiguous genes", () => {
  const tx = [
    "NM_TEST.1",
    100,
    220,
    [
      [100, 110],
      [200, 220],
    ],
    103,
    215,
    "cmpl",
    "cmpl",
  ];
  const catalog = [["TEST", "chr1", "-", 100, 220, [], [tx]]];
  const r = parseVariants(
    "NM_TEST.1:c.1T>C\nNM_TEST.1:c.15_16del\nNM_TEST.1:c.1dup\nNM_TEST.1:c.1_2insAC",
    { assembly: "hg38", catalog },
  );
  assert.equal(r.variants[0].pos, 214);
  assert.equal(r.variants[0].ref, "A");
  assert.equal(r.variants[0].alt, "G");
  assert.match(r.entries[1].error, /splice junction/);
  assert.equal(r.variants[1].copyStart, 214);
  assert.equal(r.variants[1].pos, 214);
  assert.equal(r.variants[2].alt, "GT");
  assert.equal(r.variants[2].pos, 214);
});
test("accession assembly mismatches, symbolic alleles and invalid coordinates are explicit errors", () => {
  assert.equal(parse("NC_000017.11:g.101A>C").variants[0].chr, "chr17");
  for (const s of [
    "NC_000017.10:g.101A>C",
    "chr1 0 A C",
    "chr1 101 A <DEL>",
    "chr1:g.2_4insA",
    "chr1:g.2_4delAA",
    "chr1:g.5A>A",
  ])
    assert.ok(parse(s).entries[0].error, s);
  assert.throws(
    () => resolveContigs([v(100)], [{ refName: "1", length: 100 }]),
    /outside/,
  );
  assert.equal(
    resolveContigs([v(100)], [{ refName: "1", length: 1000 }])[0].chr,
    "1",
  );
});
test("paired ends bridge a long unsequenced insert and overlapping mates count once", () => {
  const variants = [v(120), v(700, "G", 1)],
    c = makeCollector(variants, opts);
  const reads = paired(
    "bridge",
    100,
    680,
    sequence(100, 50, [[120, "C"]]),
    sequence(680, 50, [[700, "G"]]),
  );
  reads.forEach((r) => c.add(r));
  const r = summarize(variants, c.fragments, c.stats);
  assert.equal(r.fragments, 1);
  assert.equal(r.pairs[0].both, 1);
  assert.equal(r.sites[0].alt, 1);
  const overlap = makeCollector([v(120), v(130, "G", 1)], opts);
  paired(
    "overlap",
    100,
    110,
    sequence(100, 50, [
      [120, "C"],
      [130, "G"],
    ]),
    sequence(110, 50, [
      [120, "C"],
      [130, "G"],
    ]),
  ).forEach((r) => overlap.add(r));
  assert.equal(
    summarize([v(120), v(130, "G", 1)], overlap.fragments).pairs[0].both,
    1,
  );
});
test("joint frequency denominators exclude neither, OTHER and conflicting mates", () => {
  const p = {
      both: 10,
      leftOnly: 6,
      rightOnly: 4,
      neither: 8,
      other: 3,
      conflict: 2,
    },
    m = pairMetrics(p);
  assert.equal(m.together, 0.5);
  assert.equal(m.separate, 0.5);
  assert.equal(m.valid, 28);
  assert.equal(m.bGivenA, 10 / 16);
  assert.equal(
    pairMetrics({ both: 0, leftOnly: 0, rightOnly: 0, neither: 10 }).together,
    null,
  );
  const variants = [v(120), v(130, "G", 1)],
    c = makeCollector(variants, opts);
  paired(
    "conflict",
    100,
    110,
    sequence(100, 50, [
      [120, "C"],
      [130, "G"],
    ]),
    sequence(110, 50, [[130, "G"]]),
  ).forEach((r) => c.add(r));
  const result = summarize(variants, c.fragments);
  assert.equal(result.pairs[0].conflict, 1);
  assert.equal(result.pairs[0].both, 0);
});
test("gaps, low quality, missing SEQ, flags, MAPQ 255 and different samples never vote as REF", () => {
  const variant = v(120);
  for (const fields of [
    { flags: 1024 },
    { flags: 256 },
    { flags: 2048 },
    { flags: 512 },
    { mq: 255 },
    { seq: "*" },
    { seq: "A".repeat(99) },
    { sample: "other" },
  ])
    assert.equal(prepareRead(read(fields), opts), null);
  assert.equal(
    callVariant(prepareRead(read({ cigar: "20M20N80M" }), opts), variant),
    "?",
  );
  const quality = Array(100).fill(35);
  quality[20] = 255;
  assert.equal(
    callVariant(prepareRead(read({ qual: quality }), opts), variant),
    "?",
  );
  assert.equal(callVariant(prepareRead(read(), opts), v(500)), "?");
});
test("CIGAR insertion/deletion/delins, clips and sequence-free duplication are callable with flanks", () => {
  const ins = parse("chr1:g.120_121insCG").variants[0];
  let r = read({
    seq: "TTTT" + "A".repeat(20) + "CG" + "A".repeat(80),
    qual: Array(106).fill(35),
    cigar: "4S20M2I80M",
  });
  assert.equal(callVariant(prepareRead(r, opts), ins), "A");
  assert.equal(callVariant(prepareRead(read(), opts), ins), "R");
  const del = parse("chr1:g.121_122del").variants[0];
  r = read({
    seq: "A".repeat(98),
    qual: Array(98).fill(35),
    cigar: "20M2D78M",
  });
  assert.equal(callVariant(prepareRead(r, opts), del), "A");
  assert.equal(callVariant(prepareRead(read(), opts), del), "R");
  const di = parse("chr1:g.121_122delinsCG").variants[0];
  r = read({
    seq: "A".repeat(20) + "CG" + "A".repeat(78),
    qual: Array(100).fill(35),
    cigar: "20M2D2I78M",
  });
  assert.equal(callVariant(prepareRead(r, opts), di), "A");
  const dup = parse("chr1:g.119_120dup").variants[0];
  r = read({
    seq: "A".repeat(102),
    qual: Array(102).fill(35),
    cigar: "20M2I80M",
  });
  assert.equal(callVariant(prepareRead(r, opts), dup), "A");
  const edge = read({
    start: 120,
    seq: "CG" + "A".repeat(80),
    qual: Array(82).fill(35),
    cigar: "2I80M",
  });
  assert.equal(callVariant(prepareRead(edge, opts), ins), "?");
});
test("different contigs and read groups cannot create spurious local cis evidence", () => {
  const variants = [v(120), { ...v(700, "G", 1), chr: "chr2" }],
    c = makeCollector(variants, opts);
  const reads = paired(
    "x",
    100,
    680,
    sequence(100, 50, [[120, "C"]]),
    sequence(680, 50, [[700, "G"]]),
  );
  reads[0].mateChr = "chr2";
  reads[1].chr = "chr2";
  reads.forEach((r) => c.add(r));
  assert.equal(summarize(variants, c.fragments).pairs.length, 0);
  const groups = makeCollector([v(120), v(700, "G", 1)], opts);
  reads[1].chr = "chr1";
  reads[0].mateChr = "chr1";
  reads[1].rg = "other";
  reads.forEach((r) => groups.add(r));
  assert.equal(
    summarize([v(120), v(700, "G", 1)], groups.fragments).pairs.length,
    0,
  );
});
test("duplicate primary ends and inconsistent mate coordinates fail explicitly", () => {
  const c = makeCollector([v(120)], opts),
    pair = paired("bad", 100, 200, "A".repeat(50), "A".repeat(50));
  c.add(pair[0]);
  assert.throws(() => c.add({ ...pair[1], mate: 101 }), /mate coordinates/);
  const d = makeCollector([v(120)], opts);
  d.add(pair[0]);
  assert.throws(() => d.add(pair[0]), /Repeated primary/);
});
test("helper chain bridges requested sites via 5+ distinct molecules and exposes every edge", () => {
  const variants = [v(120), v(700, "C", 1)],
    options = { ...opts, helpers: true },
    result = analyzeReads(helperReads(), variants, options);
  assert.equal(result.helperCount, 2);
  assert.equal(result.pairs.filter((p) => p.a === 0 && p.b === 1).length, 0);
  const graph = helperGraph(result.variants, result.helperPairs, options),
    path = findHelperPath(graph, result.variants, 0, 1);
  assert.equal(path.status, "linked");
  assert.equal(path.sign, 1);
  assert.equal(path.helperIds.length, 2);
  assert.equal(path.edges.length, 3);
  assert.ok(path.edges.every((e) => e.count === 12 && e.total === 12));
  assert.ok(
    result.variants
      .filter((v) => v.helper)
      .every((v) => v.altCount === 12 && v.depth === 24),
  );
});
test("helper molecule count collapses repeated coordinates, excludes flagged duplicates and rejects low fractions", () => {
  const original = helperReads(),
    duplicates = original.map((r) => ({ ...r, name: r.name + "-PCR" })),
    variants = [v(120), v(700, "C", 1)];
  const helpers = discoverHelpers([...original, ...duplicates], variants, opts);
  assert.ok(helpers.every((v) => v.altCount === 12));
  const flagged = original.map((r) => ({ ...r, flags: r.flags | 1024 }));
  assert.equal(
    discoverHelpers(flagged, variants, { ...opts, duplicates: true }).length,
    0,
  );
  assert.equal(
    discoverHelpers(original, variants, { ...opts, helperSupport: 13 }).length,
    0,
  );
  assert.equal(
    discoverHelpers(original, variants, { ...opts, helperFraction: 0.5 })
      .length,
    0,
  );
  const result = analyzeReads(
    [...original, ...flagged.map((r) => ({ ...r, name: r.name + "-flagged" }))],
    variants,
    { ...opts, helpers: true, duplicates: true },
  );
  assert.ok(
    result.helperPairs.every(
      (p) => p.both + p.neither <= 12 || (p.a === 2 && p.b === 3),
    ),
  );
});
test("UMI tags distinguish same-coordinate molecules; conflicting mate calls cannot inflate support", () => {
  const original = helperReads(),
    copies = original.map((r) => ({ ...r, name: r.name + "-UMI", umi: "TAG" })),
    variants = [v(120), v(700, "C", 1)];
  assert.ok(
    discoverHelpers([...original, ...copies], variants, opts).every(
      (v) => v.altCount === 24,
    ),
  );
  const reads = paired(
    "overlap",
    100,
    110,
    sequence(100, 50, [[120, "C"]]),
    "A".repeat(50),
  );
  assert.equal(
    discoverHelpers(reads, [], { ...opts, helperSupport: 1 }).length,
    0,
  );
});
test("contradictory helper cycles are exposed and weak/monomorphic links are excluded", () => {
  const variants = [
    v(10),
    v(20, "C", 1),
    { ...v(30, "G", 2), helper: true },
    { ...v(40, "T", 3), helper: true },
  ];
  const edge = (a, b, sign = 1) => ({
    a,
    b,
    both: sign === 1 ? 5 : 0,
    neither: sign === 1 ? 5 : 0,
    leftOnly: sign === -1 ? 5 : 0,
    rightOnly: sign === -1 ? 5 : 0,
    other: 0,
    conflict: 0,
  });
  let graph = helperGraph(
    variants,
    [edge(0, 2), edge(1, 2), edge(0, 3), edge(1, 3, -1)],
    opts,
  );
  assert.equal(findHelperPath(graph, variants, 0, 1).status, "conflicting");
  graph = helperGraph(
    variants,
    [
      { ...edge(0, 2), neither: 0 },
      { ...edge(1, 2), conflict: 1 },
    ],
    opts,
  );
  assert.equal(graph.size, 0);
});
test("indexed queries stay around targets and expand only within explicit helper scope", () => {
  const variants = [v(100), v(100000, "C", 1)],
    contigs = [{ refName: "chr1", length: 200000 }];
  assert.equal(windowsForVariants(variants, contigs).length, 2);
  assert.equal(helperWindows(variants, contigs, 50000).length, 2);
  assert.equal(helperWindows(variants, contigs, 100000).length, 1);
  assert.throws(
    () =>
      helperWindows(
        [v(100), v(1999999, "C", 1), v(2999999, "C", 2)],
        [{ refName: "chr1", length: 4000000 }],
        2000000,
      ),
    /two million/,
  );
});
test("standalone specific release enforces offline CSP, embeds JS and CSS, and contains no persistence APIs", async () => {
  const html = await readFile("specific/index.html", "utf8");
  for (const tag of ["script", "style"]) {
    const s = html.match(new RegExp(`<${tag}>([\\s\\S]*?)<\\/${tag}>`))[1];
    assert.ok(
      html.includes(
        `'sha256-${createHash("sha256").update(s).digest("base64")}'`,
      ),
    );
  }
  assert.ok(html.includes("connect-src 'none'"));
  assert.doesNotMatch(
    html,
    /WebAssembly\s*\.|instantiateStreaming|localStorage|sessionStorage|indexedDB|serviceWorker|<script[^>]+src=/,
  );
});

test("wrapped transcript and ANNOVAR-style coding annotations resolve without using protein text", () => {
  const catalog = [
    [
      "TEST",
      "chr1",
      "+",
      100,
      220,
      [],
      [["NM_TEST.1", 100, 220, [[100, 220]], 103, 215, "cmpl", "cmpl"]],
    ],
  ];
  // Use a numeric synthetic accession to exercise standard accession extraction.
  catalog[0][6][0][0] = "NM_123456.1";
  const r = parseVariants(
    "NM_123456.1(TEST):c.1A>C (p.K1Q)\nTEST:NM_123456.1:exon1:c.A1C:p.K1Q",
    { assembly: "hg38", catalog },
  );
  assert.equal(r.entries.filter((e) => e.error).length, 0);
  assert.equal(r.variants.length, 1);
  assert.equal(r.entries[1].duplicate, true);
});
test("a duplicate flag on either end excludes the entire physical template in both orders", () => {
  for (const reverse of [false, true]) {
    const c = makeCollector([v(120), v(700, "G", 1)], opts),
      pair = paired(
        "one-flagged",
        100,
        680,
        sequence(100, 50, [[120, "C"]]),
        sequence(680, 50, [[700, "G"]]),
      );
    pair[1].flags |= 1024;
    if (reverse) pair.reverse();
    pair.forEach((r) => c.add(r));
    assert.equal(c.fragments.size, 0);
    assert.equal(
      discoverHelpers(pair, [], { ...opts, helperSupport: 1 }).length,
      0,
    );
  }
});
test("exact helper fraction threshold and explicit pair/decompression-discovery budgets are enforced", () => {
  const variants = [v(120), v(700, "C", 1)],
    reads = helperReads();
  assert.equal(
    discoverHelpers(reads, variants, { ...opts, helperFraction: 0.5 }).length,
    0,
  );
  assert.throws(
    () => discoverHelpers(reads, variants, opts, { maxPositions: 1 }),
    /million covered bases/,
  );
  const c = makeCollector([v(120), v(130, "G", 1)], opts);
  c.add(read());
  assert.throws(
    () => summarize([v(120), v(130, "G", 1)], c.fragments, {}, 0),
    /two million observed/,
  );
});
