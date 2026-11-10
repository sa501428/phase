import { buildOverview } from "./overview.js";
import { structuralEvidence } from "./structural.js";
import { completePairs } from "./mates.js";
import { buildFamilies } from "./families.js";
import { BamFile } from "@gmod/bam";
import { LocalFile, referenceReader, resolveReference } from "./io.js";
import { synthetic } from "./demo.js";
import {
  cropRead,
  cigarOps,
  phase,
  parseVCF,
  eventPath,
  align,
} from "./analysis.js";
let bam, ref, dict, budget, files, header;
const send = (type, data) => postMessage({ type, ...data });
self.onmessage = async ({ data: m }) => {
  try {
    if (m.type === "demo") {
      const d = synthetic();
      send("result", {
        ...d,
        result: (() => {
          const r = phase(d.reads, d.reference, d.region.start, m.opts);
          r.structural = structuralEvidence(r.reads, d.region, m.opts);
          r.families = buildFamilies(r, d.reference, d.region.start, m.opts);
          r.settings = m.opts;
          r.overview = buildOverview(r, d.region, m.opts);
          return r;
        })(),
        request: m.request,
      });
    } else if (m.type === "open") {
      files = m.files;
      budget = { bytes: 0 };
      const index = new LocalFile(files.index);
      const magic = await index.read(4);
      const kind = String.fromCharCode(...magic);
      if (
        kind !== "BAI\x01" &&
        kind !== "CSI\x01" &&
        !(magic[0] === 31 && magic[1] === 139)
      )
        throw Error("Invalid BAI/CSI index signature");
      bam = new BamFile({
        bamFilehandle: new LocalFile(files.bam, budget),
        ...(/\.csi$/i.test(files.index.name)
          ? { csiFilehandle: index }
          : { baiFilehandle: index }),
        maxCacheBytes: 16 * 1024 * 1024,
      });
      header = await bam.getHeaderText();
      if (/SO:(?!coordinate)[^\t\n]+/.test(header))
        throw Error("BAM must be coordinate sorted");
      dict = bam.indexToChr;
      ref = await referenceReader(files.reference, files.fai);
      const samples = [
        ...new Set(
          header
            .split("\n")
            .filter((x) => x.startsWith("@RG"))
            .map((x) => x.match(/\tSM:([^\t]+)/)?.[1])
            .filter(Boolean),
        ),
      ];
      send("opened", { contigs: dict, samples });
    } else if (m.type === "analyze") {
      budget.bytes = 0;
      ref.reset();
      const { region, opts } = m;
      if (
        !Number.isSafeInteger(region.start) ||
        !Number.isSafeInteger(region.end) ||
        region.start < 0 ||
        region.end <= region.start ||
        region.end - region.start > 200000
      )
        throw Error("Choose an interval between 1 and 200,000 bases");
      const d = dict.find((x) => x.refName === region.chr);
      if (!d || region.end > d.length)
        throw Error("Interval outside BAM contig");
      send("progress", {
        message: "Reading local reference and indexed alignments…",
      });
      const name = await resolveReference(ref, region.chr, d.length),
        reference = await ref.sequence(name, region.start, region.end);
      if (reference.length !== region.end - region.start)
        throw Error("Incomplete reference");
      let records = await bam.getRecordsForRange(
        region.chr,
        region.start,
        region.end,
      );
      send("progress", {
        message: "Retrieving mapped mate alignments from the local index…",
      });
      records = await completePairs(
        bam,
        records,
        region.chr,
        dict.findIndex((x) => x.refName === region.chr),
      );
      if (records.length > 100000)
        throw Error("Over 100,000 alignment records; narrow the interval");
      if (records.reduce((n, r) => n + r.seq_length, 0) > 2000000)
        throw Error(
          "Over two million sequenced bases; narrow the interval to stay within evidence memory limits",
        );
      const rgSample = new Map(
        header
          .split("\n")
          .filter((x) => x.startsWith("@RG"))
          .map((x) => [
            x.match(/\tID:([^\t]+)/)?.[1],
            x.match(/\tSM:([^\t]+)/)?.[1] || "",
          ]),
      );
      const reads = records.map((r, id) => {
        if (r.seq_length > 10000)
          throw Error("Read exceeds short-read limit (10 kb)");
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
          mateRef: r.next_refid,
          chr: dict[r.ref_id]?.refName || region.chr,
          mateChr: dict[r.next_refid]?.refName,
          sa: r.getTag("SA"),
          mateCigar: r.getTag("MC"),
          mateMq: r.getTag("MQ"),
        };
      });
      send("progress", {
        message: "Constructing fragment evidence and phase blocks…",
      });
      let supplied = [];
      if (files.vcf) {
        if (files.vcf.size > 16e6)
          throw Error("VCF exceeds 16 MiB; subset externally");
        supplied = parseVCF(
          await files.vcf.text(),
          region.chr,
          reference,
          region.start,
        );
      }
      if (m.mask?.length) {
        for (const r of reads) {
          let p = r.start,
            q = 0;
          for (const [n, op] of cigarOps(r.cigar)) {
            if ("M=X".includes(op)) {
              for (let j = 0; j < n; j++)
                if (m.mask.some(([s, e]) => p + j >= s && p + j < e))
                  r.qual[q + j] = -1;
              p += n;
              q += n;
            } else if (op === "I" || op === "S") q += n;
            else if (op === "D" || op === "N") p += n;
          }
        }
      }
      const result = phase(reads, reference, region.start, opts, supplied);
      result.structural = structuralEvidence(result.reads, region, opts);
      result.families = buildFamilies(result, reference, region.start, opts);
      result.settings = opts;
      result.overview = buildOverview(result, region, opts);
      send("result", { result, reference, region, request: m.request });
    } else if (m.type === "event") {
      const path = eventPath(m.event, m.reference, m.origin);
      const supporting = m.reads.filter(
        (r) => r.start < path.end && r.end > path.start,
      );
      const realigned = [];
      let cells = 0;
      for (const r of supporting) {
        const start = Math.max(0, r.start - path.start - 20),
          end = Math.min(
            path.alt.length,
            r.end - path.start + 20 + m.event.alt.length,
          ),
          candidate = path.alt.slice(start, end);
        if (!candidate) continue;
        const read = cropRead(r, path.start, path.end);
        if (!read) continue;
        cells += (candidate.length + 1) * (read.length + 1);
        if (cells > 20000000)
          throw Error(
            "Event realignment exceeds the total cell budget; narrow the interval. No reads were downsampled.",
          );
        const a = align(candidate, read);
        realigned.push({ id: r.id, ...a });
      }
      send("event", { path, realigned });
    }
  } catch (e) {
    send("error", {
      message:
        e.message?.replace(/https?:\/\/\S+/g, "[address]") ||
        "Unable to process input",
    });
  }
};
