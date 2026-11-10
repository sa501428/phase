import { synthetic } from "../src/demo.js";
import { gzipSync } from "node:zlib";
import { mkdir, writeFile } from "node:fs/promises";
const i32 = (n) => {
    const b = Buffer.alloc(4);
    b.writeInt32LE(n);
    return b;
  },
  str = (s) => Buffer.from(s + "\0");
export function bgzf(bytes) {
  const gz = gzipSync(bytes),
    header = Buffer.from([
      31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0,
    ]),
    out = Buffer.concat([header, gz.subarray(10)]);
  out.writeUInt16LE(out.length - 1, 16);
  return out;
}
export function fixture(d = synthetic()) {
  const reference = "N".repeat(1000) + d.reference;
  const contigs = [["synthetic", 1800], ...(d.extraContigs || [])];
  const text =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    contigs.map(([n, l]) => `@SQ\tSN:${n}\tLN:${l}\n`).join("") +
    "@RG\tID:synthetic\tSM:synthetic\n";
  const h = Buffer.concat([
    Buffer.from("BAM\x01"),
    i32(Buffer.byteLength(text)),
    Buffer.from(text),
    i32(contigs.length),
    ...contigs.flatMap(([n, l]) => [i32(n.length + 1), str(n), i32(l)]),
  ]);
  let pos = h.length;
  const records = [];
  for (const r of [...d.reads].sort((a, b) => a.start - b.start)) {
    const name = str(r.name),
      ops = [...r.cigar.matchAll(/(\d+)([MIDNSHP=X])/g)],
      cigar = Buffer.concat(
        ops.map((x) => i32(+x[1] * 16 + "MIDNSHP=X".indexOf(x[2]))),
      ),
      seq = Buffer.alloc(Math.ceil(r.seq.length / 2));
    for (let j = 0; j < r.seq.length; j++)
      seq[j >> 1] |=
        ("=ACMGRSVTWYHKDBN".indexOf(r.seq[j]) & 15) << (j % 2 ? 0 : 4);
    const fixed = Buffer.alloc(32);
    fixed.writeInt32LE(0, 0);
    fixed.writeInt32LE(r.start, 4);
    fixed.writeUInt32LE(4681 * 65536 + 60 * 256 + name.length, 8);
    fixed.writeUInt32LE(r.flags * 65536 + ops.length, 12);
    fixed.writeInt32LE(r.seq.length, 16);
    fixed.writeInt32LE(r.mate >= 0 ? (r.mateRef ?? 0) : -1, 20);
    fixed.writeInt32LE(r.mate ?? -1, 24);
    const body = Buffer.concat([
      fixed,
      name,
      cigar,
      seq,
      Buffer.from(r.qual),
      Buffer.from("RGZsynthetic\0"),
      ...(r.sa ? [Buffer.from("SAZ" + r.sa + "\0")] : []),
      ...(r.mateMq !== undefined ? [Buffer.from("MQi"), i32(r.mateMq)] : []),
    ]);
    records.push(Buffer.concat([i32(body.length), body]));
    pos += body.length + 4;
  }
  const raw = Buffer.concat([h, ...records]);
  if (raw.length > 64000)
    throw Error("Fixture too large for single BGZF block");
  const v = (n) => {
    const b = Buffer.alloc(8);
    b.writeBigUInt64LE(BigInt(n));
    return b;
  };
  const cut = Math.floor(raw.length / 2),
    block1 = bgzf(raw.subarray(0, cut)),
    block2 = bgzf(raw.subarray(cut)),
    last = block1.length * 65536 + raw.length - cut;
  const bai = Buffer.concat([
    Buffer.from("BAI\x01"),
    i32(contigs.length),
    i32(1),
    i32(4681),
    i32(1),
    v(h.length),
    v(last),
    i32(1),
    v(h.length),
    ...contigs.slice(1).flatMap(() => [i32(0), i32(0)]),
  ]);
  const csi = bgzf(
    Buffer.concat([
      Buffer.from("CSI\x01"),
      i32(14),
      i32(5),
      i32(0),
      i32(contigs.length),
      i32(1),
      i32(4681),
      v(h.length),
      i32(1),
      v(h.length),
      v(last),
      ...contigs.slice(1).map(() => i32(0)),
    ]),
  );
  const be = (n) => {
    const b = Buffer.alloc(4);
    b.writeUInt32BE(n);
    return b;
  };
  const dna = Buffer.alloc(450);
  for (let j = 0; j < reference.length; j++)
    dna[j >> 2] |=
      Math.max(0, "TCAG".indexOf(reference[j])) << ((3 - (j % 4)) * 2);
  const twobit = Buffer.concat([
    be(0x1a412743),
    be(0),
    be(1),
    be(0),
    Buffer.from([9]),
    Buffer.from("synthetic"),
    be(30),
    be(1800),
    be(1),
    be(0),
    be(1000),
    be(0),
    be(0),
    dna,
  ]);
  return {
    bam: Buffer.concat([block1, block2, bgzf(Buffer.alloc(0))]),
    bai,
    csi,
    twobit,
    fasta: Buffer.from(">synthetic\n" + reference + "\n"),
    fai: Buffer.from("synthetic\t1800\t11\t1800\t1801\n"),
  };
}
export function pairedFixture() {
  const d = synthetic();
  d.reads = [];
  for (let i = 0; i < 12; i++) {
    const a = 1050 + i,
      b = 1400 + i;
    for (const [start, mate, flags] of [
      [a, b, 97],
      [b, a, 145],
    ])
      d.reads.push({
        id: d.reads.length,
        name: "paired-" + i,
        start,
        end: start + 50,
        seq: d.reference.slice(start - 1000, start - 950),
        cigar: "50M",
        qual: Array(50).fill(35),
        flags,
        mq: 60,
        rg: "synthetic",
        sample: "synthetic",
        mate,
      });
  }
  return fixture(d);
}
export function structuralFixture() {
  const d = synthetic();
  d.extraContigs = [["chr8", 2000000]];
  d.reads = [];
  for (let i = 0; i < 9; i++) {
    const start = 1100 + i,
      split = i < 6;
    d.reads.push({
      id: i,
      name: "structural-" + i,
      start,
      end: start + (split ? 50 : 100),
      seq: split
        ? d.reference.slice(start - 1000, start - 950) +
          d.reference.slice(300, 350)
        : d.reference.slice(start - 1000, start - 900),
      cigar: split ? "50M50S" : "100M",
      flags: split ? 0 : 65,
      qual: Array(100).fill(35),
      mq: 60,
      rg: "synthetic",
      sample: "synthetic",
      mate: split ? -1 : 7000,
      mateRef: split ? -1 : 1,
      mateMq: split ? undefined : 60,
      sa:
        i < 3
          ? "chr8,5001,+,50S50M,60,0;"
          : i < 6
            ? "synthetic,1501,-,50M50S,60,0;"
            : undefined,
    });
  }
  return fixture(d);
}
export function noisyFixture() {
  const d = synthetic();
  d.reads = [];
  for (let i = 0; i < 200; i++) {
    const start = 1050 + i,
      seq = d.reference.slice(start - 1000, start - 920).split("");
    seq[20] = seq[20] === "A" ? "C" : "A";
    d.reads.push({
      id: i,
      name: "noisy-" + i,
      start,
      end: start + 80,
      seq: seq.join(""),
      qual: Array(80).fill(35),
      cigar: "80M",
      flags: 0,
      mq: 60,
      rg: "synthetic",
      sample: "synthetic",
      mate: -1,
    });
  }
  return fixture(d);
}
if (process.argv[1]?.endsWith("fixtures.mjs")) {
  await mkdir("tests/fixtures", { recursive: true });
  for (const [name, data] of [
    ["synthetic", fixture()],
    ["paired", pairedFixture()],
    ["structural", structuralFixture()],
    ["noisy", noisyFixture()],
  ])
    for (const [ext, b] of Object.entries(data))
      await writeFile(
        "tests/fixtures/" +
          name +
          "." +
          ({ fasta: "fa", fai: "fa.fai", twobit: "2bit" }[ext] || ext),
        b,
      );
  console.log(
    "Generated synthetic variant and paired-end BAM/BAI/CSI and reference fixtures",
  );
}
