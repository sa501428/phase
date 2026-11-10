import { test } from "node:test";
import assert from "node:assert/strict";
import {
  snpDistance,
  clusterPatterns,
  cutHierarchy,
  buildHierarchy,
} from "../src/hierarchy.js";
import { phase } from "../src/analysis.js";
import { synthetic } from "../src/demo.js";
const sites = [
  { index: 0, ref: "A" },
  { index: 1, ref: "C" },
  { index: 2, ref: "G" },
];
test("distance uses jointly observed SNPs and distinguishes shared ALT from shared REF", () => {
  assert.deepEqual(snpDistance(["T", "C", "?"], ["T", "T", "G"], sites), {
    distance: 0.5,
    compared: 2,
    different: 1,
    sharedAlternate: 1,
    sharedReference: 0,
  });
  assert.equal(
    snpDistance(["?", "?", "?"], ["A", "C", "G"], sites).distance,
    null,
  );
  assert.equal(
    snpDistance(["A", "C", "G"], ["A", "C", "G"], sites).sharedReference,
    3,
  );
});
test("average linkage joins closer groups first and threshold gives transparent cluster counts", () => {
  const patterns = [
      { id: 1, pattern: ["A", "C", "G"] },
      { id: 2, pattern: ["T", "C", "G"] },
      { id: 3, pattern: ["T", "T", "T"] },
    ],
    tree = clusterPatterns(patterns, sites);
  assert.equal(tree.roots.length, 1);
  assert.equal(tree.roots[0].left.groupID, 3);
  assert.deepEqual(tree.roots[0].right.leaves, [0, 1]);
  assert.ok(Math.abs(tree.roots[0].height - 5 / 6) < 1e-12);
  assert.equal(cutHierarchy(tree.roots, 0).length, 3);
  assert.equal(cutHierarchy(tree.roots, 0.4).length, 2);
  assert.equal(cutHierarchy(tree.roots, 1).length, 1);
});
test("unknown SNPs produce a forest and cannot force a merge through a third group", () => {
  const tree = clusterPatterns(
    [
      { id: 1, pattern: ["A", "?", "?"] },
      { id: 2, pattern: ["?", "C", "?"] },
      { id: 3, pattern: ["A", "C", "?"] },
    ],
    sites,
  );
  assert.equal(tree.roots.length, 2);
  assert.equal(cutHierarchy(tree.roots, 1).length, 2);
});
test("indel-only groups remain separate without SNP evidence", () => {
  const tree = clusterPatterns(
    [
      { id: 1, pattern: ["-"] },
      { id: 2, pattern: ["GG"] },
    ],
    [],
  );
  assert.equal(tree.roots.length, 2);
});
test("group frequency and site VAF have different denominators and blocks stay independent", () => {
  const d = synthetic(),
    r = phase(d.reads, d.reference, 1000, {
      mapq: 20,
      baseq: 20,
      support: 3,
      duplicates: false,
      sample: "",
    }),
    h = buildHierarchy(r);
  assert.equal(h.length, 2);
  assert.equal(h[0].groups.length, 3);
  assert.equal(h[1].groups.length, 1);
  assert.equal(h[0].informative, 36);
  assert.ok(
    h[0].groups.every(
      (g) =>
        g.fragments === 12 && g.records === 12 && g.assignedShare === 1 / 3,
    ),
  );
  const shared = h[0].groups
    .flatMap((g) => g.alleleFractions)
    .find((v) => v.pos === 1120);
  assert.equal(shared.alt, 24);
  assert.equal(shared.total, 36);
  assert.equal(shared.fraction, 2 / 3);
  assert.equal(h[0].sites.length, 2);
  assert.equal(h[1].sites.length, 0);
  assert.equal(cutHierarchy(h[0].roots, 0).length, 2);
});
test("multiallelic SNV coordinates are compared once, not counted twice", () => {
  const r = {
    variants: [
      { pos: 10, ref: "A", alt: "C" },
      { pos: 10, ref: "A", alt: "G" },
    ],
    vectors: [["C", "C"]],
    fragments: [{ reads: [0] }],
    blocks: [
      {
        id: 1,
        start: 10,
        end: 11,
        indices: [0, 1],
        groups: [{ id: 1, pattern: ["C", "C"], members: [0], compatible: 0 }],
        ambiguous: 0,
        conflicting: 0,
      },
    ],
  };
  const b = buildHierarchy(r)[0];
  assert.equal(b.sites.length, 1);
  assert.equal(b.groups[0].alleleFractions.length, 1);
});
