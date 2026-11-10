import { test } from "node:test";
import assert from "node:assert/strict";
import {
  codingPosition,
  resolveVariantLines,
  requestedPairs,
  findTranscript,
} from "../src/variant-input.js";
import { phase } from "../src/analysis.js";
import { synthetic } from "../src/demo.js";
const plus = {
    chr: "chr1",
    strand: "+",
    transcript: [
      "NM_TEST.1",
      100,
      220,
      [
        [100, 110],
        [200, 220],
      ],
      103,
      215,
    ],
  },
  minus = { ...plus, strand: "-" };
const region = { chr: "chr1", start: 100, end: 230 },
  reference = "A".repeat(130),
  resolve = (text, selected = plus) =>
    resolveVariantLines(text, { region, reference, selected });
test("compact genomic, HGVS genomic and whitespace allele tables resolve identically", () => {
  const entries = resolve("chr1:105A>G\n1:g.105A>G\nchr1 105 A G");
  assert.ok(entries.every((e) => !e.error));
  assert.equal(new Set(entries.map((e) => e.variant.id)).size, 1);
  assert.equal(entries[0].pos, 104);
});
test("coding positions use CDS start and exon-spliced coordinates", () => {
  assert.equal(codingPosition("1", plus), 103);
  assert.equal(codingPosition("7", plus), 109);
  assert.equal(codingPosition("8", plus), 200);
  assert.equal(codingPosition("-1", plus), 102);
  assert.equal(codingPosition("*1", plus), 215);
  assert.equal(codingPosition("7+2", plus), 111);
  assert.equal(codingPosition("8-2", plus), 198);
  assert.throws(() => codingPosition("2+2", plus), /boundary/);
  assert.throws(() => codingPosition("0", plus));
  assert.throws(() => codingPosition("*0", plus));
});
test("negative-strand coding substitutions reverse-complement both alleles", () => {
  assert.equal(codingPosition("1", minus), 214);
  assert.equal(codingPosition("15", minus), 200);
  assert.equal(codingPosition("16", minus), 109);
  assert.equal(codingPosition("15+2", minus), 198);
  assert.equal(codingPosition("16-2", minus), 111);
  const e = resolve("c.1T>C", minus)[0];
  assert.equal(e.pos, 214);
  assert.equal(e.ref, "A");
  assert.equal(e.alt, "G");
  assert.ok(!e.error);
});
test("simple deletion, insertion, duplication and delins notation are supported", () => {
  const entries = resolve(
    "chr1:g.105_106delAA\nchr1:g.105_106insGT\nc.2dupA\nc.2_3delinsCG",
  );
  assert.ok(
    entries.every((e) => !e.error),
    JSON.stringify(entries),
  );
  assert.equal(entries[0].ref, "AA");
  assert.equal(entries[1].pos, 105);
  assert.equal(entries[1].alt, "GT");
  assert.equal(entries[2].pos, 105);
  assert.equal(entries[3].alt, "CG");
});
test("line errors are explicit and do not discard good entries", () => {
  const entries = resolve(
    "chr1:105A>G\nchr2:105A>G\nchr1:999A>T\nchr1:105G>T\np.R1G",
  );
  assert.ok(entries[0].variant);
  assert.match(entries[1].error, /another contig/);
  assert.match(entries[2].error, /Outside/);
  assert.match(entries[3].error, /Reference allele/);
  assert.ok(entries[4].error);
  assert.match(
    resolveVariantLines("c.1A>T", { region, reference })[0].error,
    /Select a transcript/,
  );
});
test("versioned transcript selection is exact; ambiguous gene or unversioned contexts are rejected", () => {
  const catalog = [
    [
      "GENE",
      "chr1",
      "+",
      100,
      220,
      [],
      [
        plus.transcript,
        [...plus.transcript].map((x, i) => (i === 0 ? "NM_TEST.2" : x)),
      ],
    ],
  ];
  assert.equal(findTranscript("NM_TEST.1", catalog).transcript[0], "NM_TEST.1");
  assert.throws(() => findTranscript("GENE", catalog), /Ambiguous/);
  assert.throws(() => findTranscript("NM_TEST", catalog), /Ambiguous/);
});
test("direct cis observations and uncovered pairs have explicit evidence counts", () => {
  const d = synthetic(),
    entries = resolveVariantLines(
      "synthetic:1121C>A\nsynthetic:1161G>C\nsynthetic:g.1551_1558del",
      { region: d.region, reference: d.reference },
    );
  assert.ok(
    entries.every((e) => !e.error),
    JSON.stringify(entries),
  );
  const r = phase(
      d.reads,
      d.reference,
      1000,
      { baseq: 20, mapq: 20, support: 3, duplicates: false, sample: "" },
      entries.map((e) => e.variant),
    ),
    pairs = requestedPairs(entries, r);
  assert.equal(pairs[0].status, "Cis supported");
  assert.equal(pairs[0].both, 24);
  assert.equal(pairs[1].status, "Unphased");
  assert.equal(pairs[1].informative, 0);
});
test("separate patterns require an explicit diploid assumption for trans wording", () => {
  const entries = [
      { line: 1, variant: { id: "a", ref: "A", alt: "C" } },
      { line: 2, variant: { id: "b", ref: "G", alt: "T" } },
    ],
    r = {
      variants: entries.map((e) => e.variant),
      vectors: [...Array(3).fill(["C", "G"]), ...Array(3).fill(["A", "T"])],
    };
  assert.match(requestedPairs(entries, r)[0].status, /trans unresolved/);
  assert.match(
    requestedPairs(entries, r, { diploid: true })[0].status,
    /Trans supported/,
  );
  r.vectors.push(["C", "T"]);
  assert.equal(
    requestedPairs(entries, r, { diploid: true })[0].status,
    "Mixed patterns",
  );
});

test("a cDNA deletion across exons does not invent a contiguous genomic deletion", () => {
  const e = resolve("c.7_8del")[0];
  assert.match(e.error, /splice junction/);
});

test("incomplete coding annotations cannot silently define c.1", () => {
  assert.throws(
    () =>
      codingPosition("1", {
        ...plus,
        transcript: [...plus.transcript, "incmpl", "cmpl"],
      }),
    /incomplete/,
  );
});
