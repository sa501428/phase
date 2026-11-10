import { test } from "node:test";
import assert from "node:assert/strict";
import {
  phase,
  project,
  cigarOps,
  align,
  normalize,
  eventPath,
} from "../src/analysis.js";
import { synthetic } from "../src/demo.js";
import { fixture } from "../scripts/fixtures.mjs";
import { unzip, unzipChunkSegments } from "../src/bgzf.js";
import { referenceReader, LocalFile } from "../src/io.js";
const opts = { mapq: 20, baseq: 20, support: 3, duplicates: false, sample: "" };
test("synthetic SNVs, tandem insertion and deletion produce separate blocks", () => {
  const d = synthetic(),
    r = phase(d.reads, d.reference, d.region.start, opts);
  assert.equal(r.fragmentCount, 44);
  assert.equal(r.variants.length, 4);
  assert.equal(r.blocks.length, 2);
  assert.equal(r.blocks[0].groups.length, 3);
  assert.ok(r.events.some((x) => x.type === "Candidate tandem duplication"));
  assert.ok(r.events.some((x) => x.type === "Deletion"));
  assert.ok(r.blocks[0].groups.every((x) => x.count === 12));
});
test("unknown coverage never becomes reference", () => {
  const d = synthetic(),
    r = phase(d.reads, d.reference, d.region.start, opts);
  assert.ok(r.vectors.slice(0, 36).every((v) => v.at(-1) === "?"));
});
test("duplicate flags change evidence only when requested", () => {
  const d = synthetic();
  d.reads.forEach((r) => (r.flags |= 1024));
  assert.equal(phase(d.reads, d.reference, 1000, opts).fragmentCount, 0);
  assert.equal(
    phase(d.reads, d.reference, 1000, { ...opts, duplicates: true })
      .fragmentCount,
    44,
  );
});
test("overlapping mates counted once and conflicting bases become unknown", () => {
  const d = synthetic(),
    r = d.reads[0],
    a = { ...r, id: 0, flags: 65 },
    b = { ...r, id: 1, flags: 129 };
  let result = phase([a, b], d.reference, 1000, opts);
  assert.equal(result.fragmentCount, 1);
  assert.equal(result.coverage[60], 1);
  const seq = b.seq.split("");
  seq[10] = seq[10] === "A" ? "C" : "A";
  b.seq = seq.join("");
  result = phase([a, b], d.reference, 1000, opts);
  assert.equal(result.coverage[r.start + 10 - 1000], 0);
});
test("CIGAR validation, indel normalization and bounded affine alignment", () => {
  assert.throws(() => cigarOps("2Mbad"));
  assert.throws(() => align("A".repeat(2000), "A".repeat(2000)));
  assert.deepEqual(normalize(3, "", "A", "AAAAA", 0), {
    pos: 0,
    ref: "",
    alt: "A",
    id: "0:>A",
  });
  const a = align("ACGT", "ACGGT");
  assert.equal(a.a.replaceAll("-", ""), "ACGT");
  assert.equal(a.b.replaceAll("-", ""), "ACGGT");
  assert.ok(a.a.includes("-"));
  const p = eventPath({ pos: 2, ref: "", alt: "GG" }, "ACTA", 0);
  assert.equal(p.alt, "ACGGTA");
  assert.equal(p.map[2].anchor, 2);
});
test("BGZF decode and trimmed indexed chunk", async () => {
  const f = fixture(),
    raw = await unzip(f.bam);
  assert.equal(new TextDecoder().decode(raw.slice(0, 4)), "BAM\x01");
  const r = await unzipChunkSegments(f.bam, {
    minv: { blockPosition: 0, dataPosition: 4 },
    maxv: { blockPosition: 0, dataPosition: raw.length - 1 },
  });
  assert.ok(r.segments[0].length > 0);
  const bad = f.bam.slice();
  bad[25] ^= 255;
  await assert.rejects(() => unzip(bad));
});
test("FASTA range reads use FAI and reject mismatched coordinates", async () => {
  const f = fixture(),
    fa = new File([f.fasta], "x.fa"),
    fai = new File([f.fai], "x.fa.fai");
  const r = await referenceReader(fa, fai);
  assert.equal(await r.size("synthetic"), 1800);
  assert.equal(
    await r.sequence("synthetic", 1000, 1010),
    synthetic().reference.slice(0, 10),
  );
  await assert.rejects(() => r.sequence("synthetic", 0, 2000));
  await assert.rejects(() => new LocalFile(fa).read(-1));
});

test("two-bit reference decoder preserves unknown blocks and sequence", async () => {
  const f = fixture(),
    r = await referenceReader(new File([f.twobit], "x.2bit"));
  assert.equal(await r.size("synthetic"), 1800);
  assert.equal(await r.sequence("synthetic", 0, 10), "NNNNNNNNNN");
  assert.equal(
    await r.sequence("synthetic", 1000, 1020),
    synthetic().reference.slice(0, 20),
  );
});
test("anchored clip hypotheses require a unique reference anchor", async () => {
  const { clipEvent } = await import("../src/analysis.js");
  const d = synthetic(),
    p = 1300;
  const c = { pos: p, side: "right", seq: "GGG" + d.reference.slice(300, 325) };
  const e = clipEvent(c, d.reference, 1000);
  assert.ok(e);
  assert.equal(e.alt, "GGG");
  assert.equal(
    clipEvent(
      { pos: 1010, side: "right", seq: "AAAAAAAAAAAAAAAAAAAA" },
      "A".repeat(800),
      1000,
    ),
    null,
  );
});
