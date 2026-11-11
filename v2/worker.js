import { BamFile } from "@gmod/bam";
import { LocalFile, referenceReader, resolveReference } from "../src/io.js";
import { completePairs } from "../src/mates.js";
import { assemble, exportGFA } from "./assembly.js";
import { demo } from "./demo.js";
import { somaticDemo } from "./somatic-demo.js";
import { queryVariants } from "./queries.js";
let bam, ref, dict, header, budget, result;
const send = (type, data) => postMessage({ type, ...data });
self.onmessage = async ({ data: m }) => {
  try {
    if (m.type === "open") {
      result = null;
      budget = { bytes: 0 };
      ref = null;
      const index = new LocalFile(m.files.index),
        magic = await index.read(4),
        kind = String.fromCharCode(...magic);
      if (
        kind !== "BAI\x01" &&
        kind !== "CSI\x01" &&
        !(magic[0] === 31 && magic[1] === 139)
      )
        throw Error("Select a valid BAI or CSI index");
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
      if (m.files.reference)
        ref = await referenceReader(m.files.reference, m.files.fai);
      send("opened", { contigs: dict, reference: !!ref });
    } else if (m.type === "demo") {
      const d =
        m.scenario === "somatic" ? somaticDemo({ partials: true }) : demo();
      result = assemble(
        d.reads,
        d.region,
        m.reference ? d.reference : null,
        m.opts,
      );
      send("result", { result });
    } else if (m.type === "analyze") {
      if (!bam) throw Error("Open an indexed BAM first");
      const { region } = m;
      if (
        !Number.isSafeInteger(region.start) ||
        !Number.isSafeInteger(region.end) ||
        region.start < 0 ||
        region.end <= region.start ||
        region.end - region.start > 200000
      )
        throw Error(
          "Choose 1–200,000 bases (coordinates are 1-based inclusive in the form)",
        );
      const d = dict.find((x) => x.refName === region.chr);
      if (!d || region.end > d.length)
        throw Error("Interval outside the BAM contig");
      budget.bytes = 0;
      ref?.reset();
      let reference = null;
      if (ref) {
        const name = await resolveReference(ref, region.chr, d.length);
        reference = await ref.sequence(name, region.start, region.end);
        if (reference.length !== region.end - region.start)
          throw Error("Incomplete reference interval");
      }
      send("progress", {
        message: "Reading the indexed locus and retrieving same-contig mates…",
      });
      let records = await bam.getRecordsForRange(
        region.chr,
        region.start,
        region.end,
      );
      records = await completePairs(
        bam,
        records,
        region.chr,
        dict.findIndex((x) => x.refName === region.chr),
      );
      if (
        records.length > 100000 ||
        records.reduce((n, r) => n + r.seq_length, 0) > 2000000
      )
        throw Error(
          "Evidence exceeds 100,000 records or two million bases. Narrow the interval. No downsampling was performed.",
        );
      const rgSample = new Map(
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
          throw Error("A read exceeds the 10 kb short-read limit");
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
          sample: rgSample.get(rg) || "",
          mate: r.next_pos,
          chr: dict[r.ref_id]?.refName || region.chr,
          mateChr: dict[r.next_refid]?.refName,
          sa: r.getTag("SA"),
          mateCigar: r.getTag("MC"),
          mateMq: r.getTag("MQ"),
        };
      });
      send("progress", {
        message:
          "Assembling sequences, compacting the graph and threading fragments…",
      });
      result = assemble(reads, region, reference, m.opts);
      send("result", { result });
    } else if (m.type === "vaf") {
      if (!result) throw Error("Load a region first");
      result = assemble(result.reads, result.region, result.reference, {
        ...result.settings,
        minVaf: m.minVaf,
      });
      send("result", { result, preserveView: true });
    } else if (m.type === "query") {
      if (!result) throw Error("Assemble a locus first");
      send("query", {
        query: queryVariants(m.text, result, m.catalog, m.selected),
      });
    } else if (m.type === "gfa") {
      if (!result) throw Error("Assemble a locus first");
      send("gfa", { text: exportGFA(result) });
    }
  } catch (e) {
    send("error", { message: e.message || String(e) });
  }
};
