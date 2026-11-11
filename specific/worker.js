import { BamFile } from "@gmod/bam";
import { ungzip } from "pako";
import hg19 from "../references/hg19.js";
import hg38 from "../references/hg38.js";
import { LocalFile } from "../src/io.js";
import { parseVariants, resolveContigs } from "./variants.js";
import { windowsForVariants } from "./phasing.js";
import {
  analyzeReads,
  helperWindows,
  helperGraph,
  findHelperPath,
} from "./helpers.js";
let bam,
  contigs,
  samples,
  rgSamples,
  budget,
  parsed,
  fragments,
  result,
  catalog,
  catalogAssembly,
  bai,
  graph,
  loaded;
const send = (type, data = {}) => postMessage({ type, ...data });
function getCatalog(assembly) {
  if (catalogAssembly !== assembly) {
    catalog = JSON.parse(
      new TextDecoder().decode(
        ungzip(
          Uint8Array.from(atob(assembly === "hg19" ? hg19 : hg38), (c) =>
            c.charCodeAt(0),
          ),
        ),
      ),
    );
    catalogAssembly = assembly;
  }
  return catalog;
}
async function parse(text, assembly) {
  // The large transcript catalog is decompressed only for coding HGVS input.
  parsed = parseVariants(text, {
    assembly,
    catalog: /:c\./i.test(text) ? getCatalog(assembly) : [],
  });
  send("parsed", parsed);
}
async function open(files) {
  if (!files.bam || !files.index)
    throw Error("Choose both a BAM and its matching BAI (or CSI) index");
  budget = { bytes: 0 };
  const index = new LocalFile(files.index),
    magic = await index.read(4),
    kind = String.fromCharCode(...magic);
  bai = kind === "BAI\x01";
  const csi =
    kind === "CSI\x01" ||
    (magic[0] === 31 && magic[1] === 139 && /\.csi$/i.test(files.index.name));
  if (!bai && !csi)
    throw Error("Index does not have a valid BAI/CSI signature");
  bam = new BamFile({
    bamFilehandle: new LocalFile(files.bam, budget),
    ...(bai ? { baiFilehandle: index } : { csiFilehandle: index }),
    maxCacheBytes: 16 * 1024 * 1024,
  });
  const header = await bam.getHeaderText();
  if (/SO:(?!coordinate)[^\t\n]+/.test(header))
    throw Error("BAM must be coordinate sorted");
  contigs = bam.indexToChr;
  if (!contigs?.length) throw Error("BAM has no reference dictionary");
  rgSamples = new Map();
  for (const line of header.split("\n").filter((s) => s.startsWith("@RG\t"))) {
    const id = line.match(/\tID:([^\t]+)/)?.[1],
      sample = line.match(/\tSM:([^\t]+)/)?.[1];
    if (id) rgSamples.set(id, sample || "(unspecified)");
  }
  samples = [...new Set(rgSamples.values())];
  if (!samples.length) samples = ["(unspecified)"];
  send("opened", { contigs: contigs.length, samples });
}
function evidence(a, b, offset = 0) {
  const rows = [];
  let total = 0;
  for (const f of fragments?.values() || []) {
    if (!f.calls.has(a) || !f.calls.has(b)) continue;
    if (total >= offset && rows.length < 50)
      rows.push({ ...f, calls: [f.calls.get(a), f.calls.get(b)] });
    total++;
  }
  send("evidence", { a, b, offset, total, rows });
}
self.onmessage = async ({ data: m }) => {
  try {
    if (m.type === "open") {
      fragments = result = null;
      await open(m.files);
    } else if (m.type === "parse") await parse(m.text, m.assembly);
    else if (m.type === "evidence") evidence(m.a, m.b, m.offset);
    else if (m.type === "path")
      send("path", {
        a: m.a,
        b: m.b,
        path: findHelperPath(graph, result.variants, m.a, m.b),
      });
    else if (m.type === "vaf") {
      if (!loaded || !result) throw Error("Run phasing first");
      const options = { ...result.options, minVaf: m.minVaf };
      const analysis = analyzeReads(loaded.reads, loaded.variants, options);
      fragments = analysis.collector.fragments;
      const { collector, ...summary } = analysis;
      result = { ...summary, options, windows: result.windows };
      graph = helperGraph(result.variants, result.helperPairs, options);
      send("result", { result });
    } else if (m.type === "run") {
      result = fragments = null;
      if (!bam) throw Error("Open a BAM and matching index first");
      if (!parsed?.variants.length)
        throw Error("Parse at least one valid variant");
      if (parsed.entries.some((e) => e.error))
        throw Error("Correct all flagged variant lines before running");
      const options = { minVaf: 0.02, ...m.options };
      if (
        !Number.isFinite(options.minVaf) ||
        options.minVaf < 0 ||
        options.minVaf > 1
      )
        throw Error("Minimum VAF must be between 0 and 100%");
      if (!samples.includes(options.sample))
        throw Error("Choose one sample from the BAM header");
      for (const key of ["mapq", "baseq"])
        if (
          !Number.isInteger(options[key]) ||
          options[key] < 0 ||
          options[key] > 254
        )
          throw Error("Quality thresholds must be integers from 0 to 254");
      const variants = resolveContigs(parsed.variants, contigs);
      if (bai && variants.some((v) => v.pos + v.refLength + 1 >= 2 ** 29))
        throw Error(
          "A requested position is beyond BAI indexing capacity; choose a CSI index",
        );
      if (options.helpers) {
        if (
          !Number.isInteger(options.helperSupport) ||
          options.helperSupport < 5
        )
          throw Error("Helper SNPs require at least 5 unique fragments");
        if (!(options.helperFraction >= 0 && options.helperFraction < 0.5))
          throw Error(
            "Helper allele fraction must be between 0 and 50% (exclusive)",
          );
        if (
          !Number.isSafeInteger(options.helperGap) ||
          options.helperGap < 0 ||
          options.helperGap > 2000000
        )
          throw Error(
            "Helper search gap must be an integer from 0 to 2,000,000",
          );
        if (!(options.linkAgreement >= 0.9 && options.linkAgreement <= 1))
          throw Error("Helper link agreement must be between 90% and 100%");
      }
      const reads = [],
        seen = new Set(),
        windows = options.helpers
          ? helperWindows(variants, contigs, options.helperGap)
          : windowsForVariants(variants, contigs);
      budget.bytes = 0;
      let bases = 0;
      for (const [i, w] of windows.entries()) {
        send("progress", {
          message: `Reading indexed target ${i + 1} of ${windows.length} · ${seen.size.toLocaleString()} alignments`,
        });
        for (const r of await bam.getRecordsForRange(w.chr, w.start, w.end)) {
          if (seen.has(r.fileOffset)) continue;
          seen.add(r.fileOffset);
          bases += r.seq_length;
          if (seen.size > 100000 || bases > 20000000)
            throw Error(
              "More than 100,000 alignments or 20 million read bases. Analyze smaller batches; no downsampling was performed.",
            );
          if (r.seq_length > 100000)
            throw Error("An alignment exceeds the 100 kb read limit");
          const rg = String(r.getTag("RG") || ""),
            sample = rgSamples.get(rg) || "(unspecified)";
          reads.push({
            fileOffset: r.fileOffset,
            name: r.name,
            chr: contigs[r.ref_id]?.refName,
            start: r.start,
            seq: r.seq,
            qual: r.qual,
            cigar: r.CIGAR,
            flags: r.flags,
            mq: r.mq ?? 255,
            rg,
            sample,
            umi: String(r.getTag("RX") || r.getTag("UB") || ""),
            mate: r.next_pos,
            mateChr: contigs[r.next_refid]?.refName,
          });
        }
      }
      send("progress", {
        message: "Counting fragment-level allele combinations…",
      });
      loaded = { reads, variants };
      const analysis = analyzeReads(reads, variants, options);
      fragments = analysis.collector.fragments;
      const { collector, ...summary } = analysis;
      result = { ...summary, options, windows: windows.length };
      graph = helperGraph(result.variants, result.helperPairs, options);
      send("result", { result });
    }
  } catch (error) {
    send("error", { message: error.message || String(error) });
  }
};
