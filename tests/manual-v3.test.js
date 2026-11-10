import { test } from "node:test";
import assert from "node:assert/strict";
import { ManualWorkspace } from "../v3/model.js";
import {
  groupStatistics,
  manualHierarchy,
  siteCounts,
} from "../v3/statistics.js";
import { demo } from "../v2/demo.js";
import { somaticDemo } from "../v2/somatic-demo.js";
const model = (d = demo()) =>
  new ManualWorkspace({ ...d, dataset: "synthetic-data-digest" });
const select = (m, fn) =>
  (m.selection = new Set(m.fragments.filter(fn).map((f) => f.id)));
test("pairs start independently and merge only by explicit selection; overlapping mates count once", () => {
  const m = model();
  assert.equal(m.fragments.length, 52);
  assert.equal(m.rows().length, 52);
  select(m, (f) => f.name.startsWith("haplotype-1-"));
  const id = m.merge("Founder");
  assert.equal(m.groups.length, 1);
  assert.equal(m.groups[0].members.length, 12);
  assert.equal(m.rows().length, 41);
  m.expanded.add(id);
  assert.equal(m.rows().length, 53);
  const s = m.groupSummary(m.groups[0]);
  assert.equal(s.bases.find((c) => c.pos === 1160).total, 12);
  assert.equal(m.reads.length, 104);
});
test("a single contradictory SNP, including reference observations, vetoes an explicit merge atomically", () => {
  const m = model();
  select(
    m,
    (f) =>
      f.name.startsWith("haplotype-1-") || f.name.startsWith("haplotype-2-"),
  );
  const before = m.snapshot();
  assert.throws(() => m.merge("bad"), /contradictory base\/deletion/);
  assert.deepEqual(m.snapshot(), before);
  assert.equal(m.selection.size, 24);
});
test("contradictory insertions and observed absence of insertion cannot be silently merged", () => {
  const m = model();
  select(
    m,
    (f) =>
      f.name.startsWith("haplotype-2-") || f.name.startsWith("haplotype-3-"),
  );
  assert.throws(() => m.merge(), /contradictory insertion/);
});
test("a deletion and aligned bases at the same position are contradictory", () => {
  const m = model();
  select(m, (f) => f.name.startsWith("deletion-"));
  assert.throws(() => m.merge(), /contradictory base\/deletion/);
});
test("unknown gaps are not filled from the reference, and FASTA keeps disconnected blocks separate", () => {
  const d = demo(),
    a = d.reads[0],
    b = d.reads[72];
  const m = model({
    ...d,
    reads: [
      { ...a, id: 0, name: "left" },
      { ...b, id: 1, name: "right" },
    ],
  });
  m.selection = new Set([0, 1]);
  m.merge("Disconnected");
  const fasta = m.exportFasta();
  assert.equal((fasta.match(/^>/gm) || []).length, 2);
  assert.ok(!fasta.includes("N".repeat(20)));
});
test("hard/soft clips and inserted/deleted bases never shift consensus SNP columns", () => {
  const d = demo(),
    r = d.reads[24];
  r.seq = "TTTT" + r.seq + "AG";
  r.cigar = "2H4S" + r.cigar + "2S";
  r.qual = Array(r.seq.length).fill(35);
  const m = model(d);
  select(m, (f) => f.name.startsWith("haplotype-2-"));
  m.merge();
  const s = m.groupSummary(m.groups[0]);
  assert.equal(
    s.bases.find((c) => c.pos === 1120).values[0].value,
    r.seq[4 + 1120 - r.start],
  );
  assert.equal(s.clips.length, 2);
  assert.equal(s.clips[0].pos, r.start);
  assert.ok(!s.bases.some((c) => c.pos < r.start));
});
test("hide, restore, merge, ungroup and undo/redo preserve pair identity and membership", () => {
  const m = model();
  select(m, (f) => f.name.startsWith("haplotype-1-"));
  m.merge("H");
  m.selection = new Set([0, 1]);
  m.hideSelection();
  assert.equal(m.hidden.size, 2);
  assert.equal(m.groupSummary(m.groups[0]).fragments, 10);
  assert.equal(m.rows("position", null, true).length, 41);
  m.undo();
  assert.equal(m.hidden.size, 0);
  m.redo();
  assert.equal(m.hidden.size, 2);
  m.restoreAll();
  assert.equal(m.hidden.size, 0);
  m.selection = new Set([0, 1]);
  m.detach();
  assert.equal(m.groups[0].members.length, 10);
  m.undo();
  assert.equal(m.groups[0].members.length, 12);
  m.dissolve(m.groups[0].id);
  assert.equal(m.rows().length, 52);
});
test("hiding can affect visualization only when explicitly included in evidence", () => {
  const m = model();
  select(m, (f) => f.name.startsWith("haplotype-1-"));
  m.merge();
  m.selection = new Set([0]);
  m.hideSelection();
  m.setSettings({ mapq: 20, baseq: 20, excludeHidden: false });
  assert.equal(m.groupSummary(m.groups[0]).fragments, 12);
  assert.equal(
    siteCounts(
      m,
      m.fragments.map((f) => f.id),
      1120,
    ).total,
    36,
  );
});
test("VAF and group share use different, exposed denominators and deduplicate paired ends", () => {
  const m = model();
  select(m, (f) => f.name.startsWith("haplotype-2-"));
  m.merge();
  const s = groupStatistics(m, m.groups[0]),
    v = s.variants.find((v) => v.pos === 1120);
  assert.equal(s.members, 12);
  assert.equal(s.total, 52);
  assert.equal(v.altCount, 24);
  assert.equal(v.total, 36);
  assert.equal(v.localAlt, 12);
  assert.equal(v.localTotal, 12);
  m.selection = new Set([12]);
  m.hideSelection();
  const v2 = groupStatistics(m, m.groups[0]).variants.find(
    (v) => v.pos === 1120,
  );
  assert.equal(v2.total, 35);
  assert.equal(v2.altCount, 23);
});
test("assignment save/restore validates dataset, duplicate membership and imports atomically with undo", () => {
  const m = model();
  select(m, (f) => f.name.startsWith("haplotype-2-"));
  m.merge("Private");
  m.selection = new Set([0]);
  m.hideSelection();
  const saved = m.exportAssignments(),
    target = model();
  target.importAssignments(saved);
  assert.equal(target.groups[0].name, "Private");
  assert.equal(target.hidden.size, 1);
  target.undo();
  assert.equal(target.groups.length, 0);
  assert.equal(target.hidden.size, 0);
  target.redo();
  assert.equal(target.groups.length, 1);
  const bad = JSON.parse(saved);
  bad.groups[0].members.push(bad.groups[0].members[0]);
  assert.throws(
    () => target.importAssignments(JSON.stringify(bad)),
    /multiply assigned/,
  );
  bad.dataset = "other";
  assert.throws(
    () => target.importAssignments(JSON.stringify(bad)),
    /do not match/,
  );
  assert.equal(target.groups.length, 1);
});
test("ten manually constructed somatic groups retain unequal abundance and a similarity hierarchy", () => {
  const m = model(somaticDemo({ partials: true }));
  for (let i = 1; i <= 10; i++) {
    select(m, (f) => f.name.startsWith(`lineage-${i}-`));
    m.merge("Lineage " + i);
  }
  assert.equal(m.groups.length, 10);
  assert.deepEqual(
    m.groups.map((g) => g.members.length),
    [36, 24, 18, 12, 6, 30, 21, 15, 9, 3],
  );
  const t = manualHierarchy(m, 0.25);
  assert.equal(t.groups.length, 10);
  assert.ok(t.sites.length >= 8);
  assert.ok(t.clusters.length >= 2);
  assert.equal(m.rows().length, 15);
  assert.equal(m.fragments.length, 179);
});
test("hierarchy retains disconnected evidence rather than inventing matching SNP calls", () => {
  const d = demo(),
    m = model(d);
  select(m, (f) => f.name.startsWith("haplotype-2-"));
  m.merge("A");
  select(m, (f) => f.name.startsWith("deletion-1-"));
  m.merge("B");
  const t = manualHierarchy(m, 1);
  assert.equal(t.roots.length, 2);
  assert.equal(t.pairs[0][1].distance, null);
});
test("valid missing SEQ and malformed geometry stay available and contribute no fake bases", () => {
  const d = demo();
  d.reads.push(
    { ...d.reads[0], id: 104, name: "no-seq", seq: "" },
    { ...d.reads[0], id: 105, name: "bad", seq: "ACG" },
  );
  const m = model(d);
  assert.equal(m.reads.length, 106);
  assert.equal(m.fragments.length, 54);
  assert.equal(m.evidence(52).bases.size, 0);
  assert.equal(m.evidence(53).bases.size, 0);
  m.selection = new Set([52, 53]);
  m.merge("Unavailable");
  assert.equal(m.groupSummary(m.groups[0]).unknown, 2);
});
test("lowering quality thresholds can reveal conflicts without overwriting manual assignments", () => {
  const d = demo(),
    r = d.reads[0],
    q = 50;
  r.seq =
    r.seq.slice(0, q) + (r.seq[q] === "A" ? "T" : "A") + r.seq.slice(q + 1);
  r.qual[q] = 5;
  const m = model(d);
  select(m, (f) => f.name.startsWith("haplotype-1-"));
  m.merge();
  assert.equal(m.groupSummary(m.groups[0]).conflicts.length, 0);
  m.setSettings({ mapq: 20, baseq: 0, excludeHidden: true });
  assert.ok(m.groupSummary(m.groups[0]).conflicts.length);
  assert.equal(m.groups[0].members.length, 12);
  assert.throws(() => m.exportFasta(), /Resolve conflicting/);
});

test("iteratively adding selected pairs preserves a haplotype identity, while a selected subset can branch", () => {
  const m = model();
  m.selection = new Set([0, 1]);
  const id = m.merge("Founder");
  m.selection = new Set([0, 1, 2, 3]);
  assert.equal(m.merge(), id);
  assert.equal(m.groups[0].name, "Founder");
  assert.deepEqual(m.groups[0].members, [0, 1, 2, 3]);
  m.selection = new Set([0, 4]);
  const branch = m.merge("Branch");
  assert.notEqual(branch, id);
  assert.deepEqual(m.groups.find((g) => g.id === id).members, [1, 2, 3]);
  assert.deepEqual(m.groups.find((g) => g.id === branch).members, [0, 4]);
  m.undo();
  assert.equal(m.groups.length, 1);
  assert.deepEqual(m.groups[0].members, [0, 1, 2, 3]);
});
