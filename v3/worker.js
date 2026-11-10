import { BamFile } from "@gmod/bam";
import { LocalFile, referenceReader, resolveReference } from "../src/io.js";
import { completePairs } from "../src/mates.js";
import { demo } from "../v2/demo.js";
import { somaticDemo } from "../v2/somatic-demo.js";
let bam, ref, dict, header, budget;
const send = (type, data) => postMessage({ type, ...data });
async function fingerprint(d) {
  const stable = d.reads
    .map((r) => [
      r.name,
      r.sample,
      r.rg,
      r.start,
      r.end,
      r.flags,
      r.cigar,
      r.seq,
      r.mq,
      r.qual,
    ])
    .sort((a, b) => JSON.stringify(a).localeCompare(JSON.stringify(b)));
  const bytes = new TextEncoder().encode(
    JSON.stringify([d.region, d.reference, stable]),
  );
  return [...new Uint8Array(await crypto.subtle.digest("SHA-256", bytes))]
    .map((b) => b.toString(16).padStart(2, "0"))
    .join("");
}
self.onmessage = async ({ data: m }) => {
  try {
    if (m.type === "open") {
      if (!m.files.bam || !m.files.index || !m.files.reference)
        throw Error(
          "Choose BAM, BAI/CSI and a matching 2bit or indexed FASTA reference",
        );
      budget = { bytes: 0 };
      const index = new LocalFile(m.files.index),
        magic = await index.read(4),
        kind = String.fromCharCode(...magic);
      if (
        kind !== "BAI\x01" &&
        kind !== "CSI\x01" &&
        !(magic[0] === 31 && magic[1] === 139)
      )
        throw Error("Invalid BAI/CSI signature");
      bam = new BamFile({
        bamFilehandle: new LocalFile(m.files.bam, budget),
        ...(/\.csi$/i.test(m.files.index.name)
          ? { csiFilehandle: index }
          : { baiFilehandle: index }),
        maxCacheBytes: 16 * 1024 * 1024,
      });
      header = await bam.getHeaderText();
      dict = bam.indexToChr;
      if (/SO:(?!coordinate)[^\t\n]+/.test(header))
        throw Error("BAM must be coordinate sorted");
      ref = await referenceReader(m.files.reference, m.files.fai);
      send("opened", { contigs: dict });
    } else if (m.type === "demo") {
      const source =
        m.scenario === "somatic" ? somaticDemo({ partials: true }) : demo();
      const d = {
        reads: source.reads,
        reference: source.reference,
        region: source.region,
      };
      if (m.scenario !== "somatic") {
        // Legal clips on both strands: stored BAM SEQ is already oriented for
        // its alignment. The extra sequence must not shift aligned SNP columns.
        for (let i = 0; i < d.reads.length; i += 6) {
          const r = d.reads[i];
          r.seq = "TGCA" + r.seq + "CAT";
          r.cigar = "4S" + r.cigar + "3S";
          r.qual = Array(r.seq.length).fill(35);
        }
        d.reads[0].mq = 42;
        d.reads[1].mq = 58;
      }
      send("result", { data: { ...d, dataset: await fingerprint(d) } });
    } else if (m.type === "load") {
      if (!bam || !ref) throw Error("Open BAM, index and reference first");
      const { region } = m;
      if (
        !Number.isSafeInteger(region.start) ||
        !Number.isSafeInteger(region.end) ||
        region.start < 0 ||
        region.end <= region.start ||
        region.end - region.start > 200000
      )
        throw Error("Enter a gene/exon or interval of 1–200,000 bases");
      const c = dict.find((c) => c.refName === region.chr);
      if (!c || region.end > c.length)
        throw Error("Region is outside the BAM contig");
      budget.bytes = 0;
      ref.reset();
      send("progress", {
        message:
          "Reading the selected reference interval and indexed read pairs…",
      });
      const name = await resolveReference(ref, region.chr, c.length),
        reference = await ref.sequence(name, region.start, region.end);
      if (reference.length !== region.end - region.start)
        throw Error("Incomplete reference interval");
      let records = await bam.getRecordsForRange(
        region.chr,
        region.start,
        region.end,
      );
      records = await completePairs(
        bam,
        records,
        region.chr,
        dict.findIndex((c) => c.refName === region.chr),
      );
      if (
        records.length > 100000 ||
        records.reduce((n, r) => n + r.seq_length, 0) > 2000000
      )
        throw Error(
          "More than 100,000 records or two million read bases. Narrow the selected locus. No downsampling was performed.",
        );
      const samples = new Map(
        header
          .split("\n")
          .filter((s) => s.startsWith("@RG"))
          .map((s) => [
            s.match(/\tID:([^\t]+)/)?.[1],
            s.match(/\tSM:([^\t]+)/)?.[1] || "",
          ]),
      );
      const reads = records.map((r, id) => {
        if (r.seq_length > 10000)
          throw Error("Read exceeds 10 kb short-read limit");
        const rg = String(r.getTag("RG") || "");
        return {
          id,
          name: r.name,
          start: r.start,
          end: r.end,
          seq: r.seq,
          qual: [...(r.qual || [])],
          cigar: r.CIGAR,
          flags: r.flags,
          mq: r.mq ?? 255,
          rg,
          sample: samples.get(rg) || "",
          chr: dict[r.ref_id]?.refName || region.chr,
          mate: r.next_pos,
          mateChr: dict[r.next_refid]?.refName,
          sa: r.getTag("SA"),
        };
      });
      const d = { reads, reference, region };
      send("result", { data: { ...d, dataset: await fingerprint(d) } });
    }
  } catch (e) {
    send("error", { message: e.message || String(e) });
  }
};
