import { test } from "node:test";
import assert from "node:assert/strict";
import { compatibleCalls, partitionCompatible } from "../src/compatibility.js";
import { phase } from "../src/analysis.js";
import { buildFamilies } from "../src/families.js";
import { buildOverview } from "../src/overview.js";
const opts = { mapq: 20, baseq: 20, support: 3, duplicates: false, sample: "" };
const read = (id, seq, extra = {}) => ({
  id,
  name: "r" + id,
  start: 0,
  end: seq.length,
  seq,
  qual: Array(seq.length).fill(35),
  cigar: seq.length + "M",
  flags: 0,
  mq: 60,
  rg: "r",
  sample: "s",
  mate: -1,
  ...extra,
});
function model(reads, ref, settings = opts) {
  const r = phase(reads, ref, 0, settings);
  r.families = buildFamilies(r, ref, 0, opts);
  r.overview = buildOverview(r, { start: 0, end: ref.length }, opts);
  return r;
}
function verifyRows(r) {
  for (const group of r.overview) {
    const calls = new Map();
    for (const member of group.members)
      for (const id of member.readIds) {
        const read = r.reads.find((r) => r.id === id);
        for (let p = read.start; p < read.end; p++) {
          const base = read.seq[p - read.start];
          if (calls.has(p))
            assert.equal(calls.get(p), base, `clashing group at ${p}`);
          calls.set(p, base);
        }
      }
  }
}
test("one SNP contradiction is an absolute veto, reference observations included", () => {
  assert.equal(
    compatibleCalls(new Map([[10, "A"]]), new Map([[10, "C"]])),
    false,
  );
  assert.equal(
    compatibleCalls(new Map([[10, "A"]]), new Map([[20, "C"]])),
    true,
  );
});
test("partial observations cannot bridge incompatible subgroups", () => {
  const items = [
      new Map([
        [0, "A"],
        [1, "C"],
      ]),
      new Map([
        [0, "G"],
        [1, "C"],
      ]),
      new Map([[1, "C"]]),
    ],
    parts = partitionCompatible(items, (m) => m);
  assert.equal(parts.length, 3);
  for (const p of parts)
    for (const item of p.items) assert.ok(compatibleCalls(p.calls, item));
});
test("singleton SNP below discovery threshold is neither assigned nor collapsed with contrary reads", () => {
  const ref = "A".repeat(30),
    base = "A".repeat(10) + "C" + "A".repeat(19),
    reads = Array.from({ length: 6 }, (_, id) => read(id, base));
  reads.push(read(6, base.slice(0, 20) + "G" + base.slice(21)));
  const r = model(reads, ref);
  assert.equal(r.variants.length, 1);
  const g = r.blocks[0].groups[0];
  assert.equal(g.members.length, 6);
  assert.ok(!g.members.includes(6));
  assert.ok(
    !r.overview.some(
      (row) => row.members.some((m) => m.id === 6) && row.members.length > 1,
    ),
  );
  verifyRows(r);
});
test("low-quality visible base disagreements still split display collapse with suppression disabled", () => {
  const ref = "A".repeat(20),
    reads = [
      read(0, ref),
      read(1, ref),
      read(2, ref),
      read(3, "A".repeat(10) + "C" + "A".repeat(9), {
        qual: Array(20).fill(5),
      }),
    ],
    r = model(reads, ref, { ...opts, minVaf: 0 });
  assert.ok(r.overview.length > 1);
  verifyRows(r);
});
test("compatible partials and duplicate records remain collapsible", () => {
  const ref = "A".repeat(20),
    seq = "A".repeat(10) + "C" + "A".repeat(9),
    reads = [
      read(0, seq),
      read(1, seq),
      read(2, seq),
      read(3, "A".repeat(5) + "C" + "A".repeat(4), { start: 5, end: 15 }),
      read(4, seq, { flags: 1024 }),
    ],
    r = model(reads, ref);
  assert.equal(r.overview.length, 1);
  assert.equal(r.overview[0].count, 5);
  verifyRows(r);
});

test("identical discordant mate pairs stay separate in exact and overview modes", () => {
  const reads = [
      read(0, "AAAA", { name: "a", flags: 65 }),
      read(1, "ACAA", { name: "a", flags: 129 }),
      read(2, "AAAA", { name: "b", flags: 65 }),
      read(3, "ACAA", { name: "b", flags: 129 }),
    ],
    r = model(reads, "AAAA");
  assert.equal(r.families.length, 2);
  assert.equal(r.overview.length, 2);
  assert.ok(r.overview.every((g) => g.discordantPair && g.count === 1));
});
