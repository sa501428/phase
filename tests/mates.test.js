import { test } from "node:test";
import assert from "node:assert/strict";
import { completePairs } from "../src/mates.js";
const r = (name, start, flags, next_pos, id, rg = "a") => ({
  name,
  start,
  flags,
  next_pos,
  next_refid: 0,
  fileOffset: id,
  getTag: () => rg,
});
test("indexed mate retrieval keeps both ends, rejects unrelated names and read groups", async () => {
  const first = r("pair", 10, 65, 100, 1),
    mate = r("pair", 100, 129, 10, 2);
  let n = 0;
  const bam = {
    getRecordsForRange: async (chr, start, end) => {
      n++;
      assert.equal(chr, "chr1");
      assert.equal(start, 0);
      return [
        mate,
        r("other", 100, 129, 10, 3),
        r("pair", 100, 129, 10, 4, "wrong"),
      ];
    },
  };
  const records = await completePairs(bam, [first], "chr1", 0);
  assert.equal(records.length, 2);
  assert.equal(records[1], mate);
  assert.equal(n, 1);
});
test("existing mates cause no additional index requests", async () => {
  const records = [r("pair", 10, 65, 100, 1), r("pair", 100, 129, 10, 2)];
  assert.equal(
    (
      await completePairs(
        {
          getRecordsForRange: () => {
            throw Error("Unexpected lookup");
          },
        },
        records,
        "chr1",
        0,
      )
    ).length,
    2,
  );
});
test("mate limits fail explicitly rather than sampling requests", async () => {
  const records = Array.from({ length: 513 }, (_, i) =>
    r("p" + i, 1, 65, i * 1000, 1000 + i),
  );
  await assert.rejects(
    () => completePairs({}, records, "chr1", 0),
    /No reads were downsampled/,
  );
});
