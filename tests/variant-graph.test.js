import { test } from "node:test";
import assert from "node:assert/strict";
import { ManualWorkspace } from "../v3/model.js";
import {
  buildVariantGraph,
  layoutVariantGraph,
  edgeAppearance,
} from "../v3/variant-graph.js";
import { demo } from "../v2/demo.js";
import { somaticDemo } from "../v2/somatic-demo.js";
const model = (d) => new ManualWorkspace({ ...d, dataset: "graph-test" });
const index = (g, pos, kind = "SNV") =>
  g.nodes.findIndex((n) => n.pos === pos && n.kind === kind);
const edge = (g, a, b) =>
  g.edges.find((e) => e.a === Math.min(a, b) && e.b === Math.max(a, b));
test("graph nodes are alleles, including a single multi-base deletion; mates count once and disconnected blocks stay unlinked", () => {
  const m = model(demo()),
    g = buildVariantGraph(m);
  assert.equal(g.fragments, 52);
  assert.equal(g.nodes.length, 4);
  const a = index(g, 1120),
    b = index(g, 1160),
    i = g.nodes.findIndex((n) => n.kind === "insertion"),
    d = g.nodes.findIndex((n) => n.kind === "deletion");
  assert.ok(a >= 0 && b >= 0 && i >= 0 && d >= 0);
  assert.equal(g.nodes[d].ref.length, 8);
  assert.equal(g.nodes[a].altCount, 24);
  assert.equal(g.nodes[a].total, 36);
  const linked = edge(g, a, b);
  assert.equal(linked.both, 24);
  assert.equal(linked.callable, 36);
  assert.equal(linked.neither, 12);
  assert.equal(linked.unknown, 16);
  assert.equal(edge(g, a, i).both, 12);
  assert.equal(edge(g, a, i).leftOnly, 12);
  assert.equal(edge(g, a, i).neither, 12);
  assert.ok(g.edges.every((e) => e.a !== d && e.b !== d));
  assert.ok(g.nodes.every((n) => n.carriers.length === n.altCount));
});
test("somatic graph shows nested and sibling evidence without assuming diploidy or manufacturing a full clone", () => {
  const g = buildVariantGraph(model(somaticDemo()));
  const founder = index(g, 1050),
    otherFounder = index(g, 1090),
    branch = index(g, 1130),
    child = index(g, 1410);
  assert.equal(edge(g, founder, branch).both, 30);
  assert.equal(edge(g, branch, child).both, 6);
  assert.equal(edge(g, founder, otherFounder), undefined);
  assert.equal(
    index(g, 1450),
    -1,
    "private allele below 2% is not rendered as a node",
  );
  assert.equal(edge(g, branch, child).callable, 174);
  const largest = Math.max(...g.edges.map((e) => e.frequency));
  assert.ok(
    edgeAppearance(edge(g, founder, branch).frequency, largest).width >
      5 * edgeAppearance(edge(g, branch, child).frequency, largest).width,
  );
});
test("live VAF changes remove and restore graph variants and all their edges", () => {
  const m = model(somaticDemo());
  m.setSettings({ ...m.settings, minVaf: 0.01 });
  let g = buildVariantGraph(m);
  assert.ok(index(g, 1450) >= 0);
  m.setSettings({ ...m.settings, minVaf: 0.05 });
  g = buildVariantGraph(m);
  assert.equal(index(g, 1450), -1);
  assert.equal(index(g, 1410), -1);
  assert.ok(g.edges.every((e) => g.nodes[e.a] && g.nodes[e.b]));
  m.setSettings({ ...m.settings, minVaf: 0.01 });
  assert.ok(index(buildVariantGraph(m), 1450) >= 0);
});
function r(id, name, start, seq, extra = {}) {
  return {
    id,
    name,
    start,
    end: start + seq.length,
    seq,
    qual: Array(seq.length).fill(35),
    cigar: seq.length + "M",
    mq: 60,
    flags: 0,
    rg: "g",
    sample: "s",
    chr: "chr1",
    ...extra,
  };
}
test("paired ends bridge sites, while a one-sided ALT call never becomes a negative observation at its missing site", () => {
  const reads = [];
  for (let f = 0; f < 10; f++)
    for (const [p, flags] of [
      [0, 99],
      [200, 147],
    ]) {
      let seq = Array(50).fill("A");
      if (p === 0 && f < 8) seq[20] = "C";
      if (p === 200 && (f < 6 || f >= 8)) seq[20] = "G";
      reads.push(r(reads.length, "pair" + f, p, seq.join(""), { flags }));
    }
  reads.push(
    r(reads.length, "one-sided", 200, "A".repeat(20) + "G" + "A".repeat(29)),
  );
  const g = buildVariantGraph(
    model({
      reads,
      reference: "A".repeat(300),
      region: { chr: "chr1", start: 0, end: 300 },
    }),
  );
  const e = g.edges[0];
  assert.equal(e.both, 6);
  assert.equal(e.leftOnly, 2);
  assert.equal(e.rightOnly, 2);
  assert.equal(e.neither, 0);
  assert.equal(e.callable, 10);
  assert.equal(e.unknown, 1);
});
test("contradictory mates and other ALT alleles do not enter binary pair frequencies", () => {
  const reference = "A".repeat(80),
    reads = [];
  for (let i = 0; i < 40; i++) {
    let s = Array(80).fill("A");
    s[20] = i < 20 ? "C" : i < 30 ? "G" : "A";
    s[40] = i < 25 ? "T" : "A";
    reads.push(r(reads.length, "f" + i, 0, s.join(""), { flags: 99 }));
  }
  let conflict = Array(80).fill("A");
  conflict[40] = "T";
  reads.push(r(reads.length, "f0", 0, conflict.join(""), { flags: 147 }));
  const g = buildVariantGraph(
    model({ reads, reference, region: { chr: "chr1", start: 0, end: 80 } }),
  );
  const c = g.nodes.find((n) => n.pos === 20 && n.alt === "C"),
    gg = g.nodes.find((n) => n.pos === 20 && n.alt === "G"),
    t = g.nodes.find((n) => n.pos === 40);
  assert.equal(c.altCount, 19);
  assert.equal(c.conflicting, 1);
  assert.equal(c.other, 10);
  const e = edge(g, c.index, t.index);
  assert.equal(e.both, 19);
  assert.equal(e.callable, 29);
  assert.equal(e.unknown, 11);
  assert.equal(edge(g, c.index, gg.index), undefined);
});
test("equivalent repeat-shifted indels create one allele node and do not contradict their own alignment", () => {
  const reads = [
    r(0, "a", 0, "A".repeat(98), { cigar: "30M2D68M", end: 100 }),
    r(1, "b", 0, "A".repeat(98), { cigar: "40M2D58M", end: 100 }),
    r(2, "ref", 0, "A".repeat(100)),
  ];
  const g = buildVariantGraph(
    model({
      reads,
      reference: "A".repeat(100),
      region: { chr: "chr1", start: 0, end: 100 },
    }),
  );
  assert.equal(g.nodes.length, 1);
  assert.equal(g.nodes[0].kind, "deletion");
  assert.equal(g.nodes[0].altCount, 2);
  assert.equal(g.nodes[0].total, 3);
  assert.equal(g.nodes[0].conflicting, 0);
});
test("duplicate and hidden pair filters change graph support and denominators; manual merging never adds votes", () => {
  const m = model(demo());
  let g = buildVariantGraph(m),
    before = g.nodes.find((n) => n.pos === 1120).altCount;
  m.selection = new Set(
    m.fragments
      .filter((f) => f.name.startsWith("haplotype-2"))
      .map((f) => f.id),
  );
  m.merge();
  assert.equal(
    buildVariantGraph(m).nodes.find((n) => n.pos === 1120).altCount,
    before,
  );
  m.selection = new Set(m.groups[0].members);
  m.hideSelection();
  assert.equal(
    buildVariantGraph(m).nodes.find((n) => n.pos === 1120).altCount,
    12,
  );
  m.restoreAll();
  assert.equal(
    buildVariantGraph(m).nodes.find((n) => n.pos === 1120).altCount,
    before,
  );
  const dup = demo();
  for (const rr of dup.reads.filter((r) => r.name.startsWith("haplotype-2")))
    rr.flags |= 1024;
  const filtered = model(dup);
  assert.equal(
    buildVariantGraph(filtered).nodes.find((n) => n.pos === 1120).altCount,
    12,
  );
  filtered.setSettings({ ...filtered.settings, includeDuplicates: true });
  assert.equal(
    buildVariantGraph(filtered).nodes.find((n) => n.pos === 1120).altCount,
    24,
  );
});
test("large-contrast edge encoding and layout are deterministic, finite, and show isolated nodes", () => {
  const weak = edgeAppearance(0.01, 0.4),
    strong = edgeAppearance(0.4, 0.4);
  assert.ok(strong.width > 15 * weak.width);
  assert.ok(strong.opacity > weak.opacity * 3);
  const g = buildVariantGraph(model(demo())),
    a = layoutVariantGraph(g.nodes, g.edges),
    b = layoutVariantGraph(g.nodes, g.edges);
  assert.deepEqual(a, b);
  assert.equal(a.length, g.nodes.length);
  assert.ok(a.every((p) => Number.isFinite(p.x) && Number.isFinite(p.y)));
  assert.deepEqual(layoutVariantGraph([], []), []);
  assert.match(
    buildVariantGraph(model(demo()), { maxVariants: 2 }).error,
    /All reads remain available/,
  );
});
