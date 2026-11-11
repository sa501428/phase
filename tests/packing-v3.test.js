import { test } from "node:test";
import assert from "node:assert/strict";
import { ManualWorkspace } from "../v3/model.js";
import { packPileup, pileupHit, fragmentSpan } from "../v3/layout.js";
import { siteCounts } from "../v3/statistics.js";
const region = { chr: "synthetic", start: 1000, end: 1800 };
function data(specs) {
  const reads = [];
  for (const [name, start, end, flags = 0] of specs) {
    const n = end - start;
    reads.push({
      id: reads.length,
      name,
      start,
      end,
      seq: "A".repeat(n),
      cigar: n + "M",
      qual: Array(n).fill(35),
      flags,
      mq: 60,
      chr: "synthetic",
    });
  }
  return { reads, region, reference: "A".repeat(800), dataset: "packed" };
}
const layout = (m, opts = {}) =>
  packPileup(m, m.rows(), {
    view: region,
    width: 1000,
    left: 112,
    readHeight: 18,
    ...opts,
  });
test("nonoverlapping templates share one lane without merging their identities", () => {
  const m = new ManualWorkspace(
      data([
        ["a", 1010, 1080],
        ["b", 1150, 1220],
        ["c", 1350, 1420],
      ]),
    ),
    p = layout(m);
  assert.equal(p.laneCount, 1);
  assert.equal(p.placements.length, 3);
  assert.equal(new Set(p.placements.map((p) => p.y)).size, 1);
  for (const q of p.placements) {
    const hit = pileupHit(p, (q.x1 + q.x2) / 2, q.y + 9);
    assert.equal(hit.key, q.key);
  }
  assert.equal(m.groups.length, 0);
});
test("overlapping templates occupy different lanes, including their paired insert spans", () => {
  const d = data([
      ["pair", 1010, 1050, 65],
      ["pair", 1200, 1240, 129],
      ["inside-insert", 1100, 1150],
      ["later", 1300, 1380],
    ]),
    m = new ManualWorkspace(d),
    p = layout(m);
  assert.equal(p.laneCount, 2);
  assert.equal(m.fragments.length, 3);
  assert.notEqual(
    p.placements.find((p) => p.key === "f:0").y,
    p.placements.find((p) => p.key === "f:1").y,
  );
  assert.ok(
    [
      p.placements.find((p) => p.key === "f:0").y,
      p.placements.find((p) => p.key === "f:1").y,
    ].includes(p.placements.find((p) => p.key === "f:2").y),
  );
});
test("soft clip extents and pixel padding prevent visually overlapping reads sharing a lane", () => {
  const d = data([
    ["clipped", 1050, 1100],
    ["near", 1110, 1160],
  ]);
  d.reads[0].seq += "A".repeat(25);
  d.reads[0].cigar += "25S";
  d.reads[0].qual = Array(75).fill(35);
  const m = new ManualWorkspace(d);
  assert.equal(fragmentSpan(m, m.fragments[0]).end, 1125);
  assert.equal(layout(m).laneCount, 2);
});
test("packed display contains only current-view items but every pair remains retained for panning", () => {
  const m = new ManualWorkspace(
    data([
      ["a", 1010, 1080],
      ["b", 1400, 1470],
    ]),
  );
  assert.equal(
    layout(m, { view: { start: 1000, end: 1200 } }).placements.length,
    1,
  );
  assert.equal(
    layout(m, { view: { start: 1300, end: 1500 } }).placements[0].key,
    "f:1",
  );
  assert.equal(m.fragments.length, 2);
});
test("duplicate-flagged pairs are filtered from pileup, coverage and VAF by default and can be enabled", () => {
  const m = new ManualWorkspace(
    data([
      ["a", 1010, 1080],
      ["duplicate", 1010, 1080, 1024],
    ]),
  );
  assert.equal(m.fragments[1].duplicate, true);
  assert.equal(layout(m).placements.length, 1);
  assert.equal(m.coverage()[20], 1);
  assert.equal(siteCounts(m, [0, 1], 1020).total, 1);
  m.setSettings({ mapq: 20, baseq: 20, includeDuplicates: true });
  assert.equal(layout(m).placements.length, 2);
  assert.equal(m.coverage()[20], 2);
  assert.equal(siteCounts(m, [0, 1], 1020).total, 2);
  assert.equal(m.reads.length, 2);
});
test("one duplicate-flagged end filters the whole physical pair", () => {
  const m = new ManualWorkspace(
    data([
      ["pair", 1010, 1080, 65],
      ["pair", 1100, 1170, 129 | 1024],
    ]),
  );
  assert.equal(m.fragments.length, 1);
  assert.equal(m.rows().length, 0);
  assert.equal(m.activeIDs([0]).length, 0);
});
test("manual consensus tracks stay separate while expanded members pack within their own track", () => {
  const m = new ManualWorkspace(
    data([
      ["a", 1010, 1080],
      ["b", 1150, 1220],
      ["c", 1350, 1420],
    ]),
  );
  m.selection = new Set([0, 1]);
  const id = m.merge("H1");
  m.expanded.add(id);
  const p = layout(m);
  assert.equal(p.placements[0].header, true);
  const children = p.placements.filter((p) => p.key.startsWith("c:"));
  assert.equal(children.length, 2);
  assert.equal(children[0].y, children[1].y);
  assert.ok(children[0].y > p.placements[0].y);
});
test("packing assigns every displayed item once and never overlaps items in a lane at high depth", () => {
  const specs = Array.from({ length: 1500 }, (_, i) => [
      "r" + i,
      1000 + ((i * 17) % 650),
      1050 + ((i * 17) % 650),
    ]),
    m = new ManualWorkspace(data(specs)),
    p = layout(m);
  assert.equal(p.placements.length, 1500);
  assert.equal(new Set(p.placements.map((p) => p.key)).size, 1500);
  const lanes = new Map();
  for (const q of p.placements) {
    if (!lanes.has(q.y)) lanes.set(q.y, []);
    lanes.get(q.y).push(q);
  }
  for (const lane of lanes.values()) {
    lane.sort((a, b) => a.x1 - b.x1);
    for (let i = 1; i < lane.length; i++)
      assert.ok(lane[i - 1].x2 + 8 <= lane[i].x1);
  }
});
