import { cigarOps } from "./analysis.js";

function geometry(start, cigar, strand) {
  const ops = cigarOps(cigar);
  let ref = 0,
    query = 0,
    leading = 0,
    aligned = 0;
  for (const [n, op] of ops) {
    if ("MDN=X".includes(op)) ref += n;
    if ("MISH=X".includes(op)) query += n;
    if ("MI=X".includes(op)) aligned += n;
  }
  for (const [n, op] of ops) {
    if (op === "S" || op === "H") leading += n;
    else break;
  }
  if (!ref || !aligned || query > 100000)
    throw Error("Invalid split alignment geometry");
  const qStart = strand === "-" ? query - leading - aligned : leading;
  return {
    start,
    end: start + ref,
    qStart,
    qEnd: qStart + aligned,
    query,
    strand,
  };
}
export function parseSA(text) {
  const alignments = [];
  let invalid = 0;
  if (!text) return { alignments, invalid };
  if (typeof text !== "string" || text.length > 65536)
    return { alignments, invalid: 1 };
  for (const entry of text.split(";").filter(Boolean)) {
    try {
      const fields = entry.split(",");
      if (fields.length !== 6) throw Error();
      const [chr, position, strand, cigar, mq, nm] = fields,
        pos = Number(position) - 1,
        mapq = Number(mq);
      if (
        !chr ||
        !/^\d+$/.test(position) ||
        !/^\d+$/.test(mq) ||
        !/^\d+$/.test(nm) ||
        !Number.isSafeInteger(pos) ||
        pos < 0 ||
        !["+", "-"].includes(strand) ||
        !Number.isInteger(mapq) ||
        mapq < 0 ||
        mapq > 255 ||
        !Number.isInteger(+nm) ||
        +nm < 0
      )
        throw Error();
      alignments.push({
        chr,
        pos,
        strand,
        cigar,
        mapq,
        ...geometry(pos, cigar, strand),
      });
    } catch {
      invalid++;
    }
  }
  return { alignments, invalid };
}
const identity = (r) => JSON.stringify([r.sample, r.rg, r.name]);
// No partner-region queries: all annotations come from records that overlap the
// selected region, their existing supplementary segments, and mate/SA metadata.
export function structuralEvidence(reads, region, opts) {
  const observations = [],
    warnings = [],
    byRead = new Map();
  let invalidSA = 0;
  for (const r of reads) {
    const key = JSON.stringify([r.sample, r.rg, r.name, r.flags & 192]);
    if (!byRead.has(key)) byRead.set(key, []);
    byRead.get(key).push(r);
  }
  for (const r of reads) {
    const chr = r.chr || region.chr;
    if (
      chr !== region.chr ||
      r.start >= region.end ||
      r.end <= region.start ||
      r.flags & (4 | 256)
    )
      continue;
    if (opts.sample && r.sample !== opts.sample) continue;
    const localOK =
      !(r.flags & 512) &&
      !(!opts.duplicates && r.flags & 1024) &&
      r.mq !== 255 &&
      r.mq >= opts.mapq;
    const strand = r.flags & 16 ? "-" : "+";
    let source;
    try {
      source = geometry(r.start, r.cigar, strand);
    } catch {
      continue;
    }
    const parsed = parseSA(r.sa);
    invalidSA += parsed.invalid;
    const targets = [
      ...parsed.alignments.map((a) => ({ ...a, via: "SA tag" })),
    ];
    for (const other of byRead.get(
      JSON.stringify([r.sample, r.rg, r.name, r.flags & 192]),
    ) || []) {
      if (
        other.id === r.id ||
        !(r.flags & 2048 || other.flags & 2048) ||
        other.flags & (4 | 256)
      )
        continue;
      try {
        targets.push({
          chr: other.chr || region.chr,
          pos: other.start,
          cigar: other.cigar,
          mapq: other.mq,
          ...geometry(other.start, other.cigar, other.flags & 16 ? "-" : "+"),
          via: "supplementary record",
        });
      } catch {}
    }
    for (const t of targets) {
      if (t.chr === chr && t.start === source.start && t.strand === strand)
        continue;
      const overlap = Math.max(
        0,
        Math.min(source.qEnd, t.qEnd) - Math.max(source.qStart, t.qStart),
      );
      const smallerSpan = Math.min(
        source.qEnd - source.qStart,
        t.qEnd - t.qStart,
      );
      const complementary =
        source.query === t.query &&
        overlap <= Math.min(10, Math.floor(smallerSpan * 0.2)) &&
        smallerSpan - overlap >= 12;
      const sourceFirst = source.qStart < t.qStart;
      const localPos = sourceFirst
        ? strand === "+"
          ? source.end - 1
          : source.start
        : strand === "+"
          ? source.start
          : source.end - 1;
      const partnerPos = sourceFirst
        ? t.strand === "+"
          ? t.start
          : t.end - 1
        : t.strand === "+"
          ? t.end - 1
          : t.start;
      let kind, type;
      if (!complementary) {
        kind = "ambiguous-split";
        type = "Ambiguous supplementary alignment";
      } else if (t.chr !== chr) {
        kind = "translocation-split";
        type = "Possible translocation · split alignment";
      } else if (strand !== t.strand) {
        kind = "inversion-split";
        type = "Inversion-like split alignment";
      } else {
        kind = "regional-split";
        type = "Split alignment to another region";
      }
      observations.push({
        kind,
        type,
        local: {
          chr,
          pos: localPos,
          strand,
          start: source.start,
          end: source.end,
        },
        partner: {
          chr: t.chr,
          pos: partnerPos,
          strand: t.strand,
          start: t.start,
          end: t.end,
        },
        source: t.via,
        record: r.id,
        fragment: identity(r),
        eligible: localOK && t.mapq !== 255 && t.mapq >= opts.mapq,
        localOK,
        partnerMapq: t.mapq,
        duplicate: !!(r.flags & 1024),
        note: complementary
          ? "Complementary segments on the same sequenced read; breakpoint coordinates remain approximate."
          : `Query segments overlap by ${overlap} bases or have incompatible lengths. Alternative/repetitive mapping is possible; not treated as a rearrangement call.`,
      });
    }
    // Mate flags describe the other end even when its record is outside the
    // selected region or on a chromosome we have not loaded.
    if (!(r.flags & 1) || r.flags & (8 | 2048) || !r.mateChr || r.mate < 0)
      continue;
    const mateStrand = r.flags & 32 ? "-" : "+",
      cross = r.mateChr !== chr,
      distance = cross ? null : Math.abs(r.mate - r.start),
      sameStrand = strand === mateStrand;
    if (!cross && !sameStrand && distance <= 10000) continue;
    const kind = cross
      ? "translocation-pair"
      : sameStrand
        ? "inversion-pair"
        : "distant-pair";
    const type = cross
      ? "Possible translocation · mate on another chromosome"
      : sameStrand
        ? "Orientation-discordant pair · possible inversion"
        : "Mate aligned to a distant region";
    const opposite = r.flags & 64 ? 128 : r.flags & 128 ? 64 : 0;
    const loadedMate = (
      byRead.get(JSON.stringify([r.sample, r.rg, r.name, opposite])) || []
    ).find(
      (other) =>
        other.start === r.mate &&
        (other.chr || region.chr) === r.mateChr &&
        !(other.flags & (4 | 256 | 2048)),
    );
    const partnerMapq = r.mateMq ?? loadedMate?.mq;
    let partnerEnd = r.mate + 1;
    try {
      if (r.mateCigar)
        partnerEnd = geometry(r.mate, r.mateCigar, mateStrand).end;
    } catch {}
    observations.push({
      kind,
      type,
      local: { chr, pos: r.start, strand, start: r.start, end: r.end },
      partner: {
        chr: r.mateChr,
        pos: r.mate,
        strand: mateStrand,
        start: r.mate,
        end: partnerEnd,
      },
      source: "mate metadata",
      record: r.id,
      fragment: identity(r),
      eligible:
        localOK &&
        partnerMapq !== undefined &&
        partnerMapq !== 255 &&
        partnerMapq >= opts.mapq,
      localOK,
      partnerMapq,
      duplicate: !!(r.flags & 1024),
      note: cross
        ? "Interchromosomal mate link; mate metadata alone does not resolve the breakpoint."
        : sameStrand
          ? "Same-strand mate orientation. This is only an inversion clue for libraries expected to have opposite-strand pairs; no library orientation is assumed or calibrated."
          : "Mate start is more than 10 kb away. This is a display threshold, not an assay-calibrated structural variant call.",
    });
  }
  // Nearby anchors are clustered within 100 bp of the first observation; keep
  // full ranges and exact supporting records so the tolerance stays transparent.
  const buckets = new Map(),
    clusters = [];
  for (const o of observations) {
    const base = JSON.stringify([
        o.kind,
        o.local.chr,
        o.partner.chr,
        o.local.strand,
        o.partner.strand,
      ]),
      bin = Math.floor(o.local.pos / 100);
    let cluster;
    for (let b = bin - 1; b <= bin + 1 && !cluster; b++)
      for (const c of buckets.get(base + ":" + b) || [])
        if (
          Math.abs(c.anchorLocal - o.local.pos) <= 100 &&
          Math.abs(c.anchorPartner - o.partner.pos) <= 100
        ) {
          cluster = c;
          break;
        }
    if (!cluster) {
      cluster = {
        id: clusters.length + 1,
        kind: o.kind,
        type: o.type,
        local: { ...o.local },
        partner: { ...o.partner },
        anchorLocal: o.local.pos,
        anchorPartner: o.partner.pos,
        localRange: [o.local.pos, o.local.pos],
        partnerRange: [o.partner.pos, o.partner.pos],
        observations: [],
      };
      clusters.push(cluster);
      const key = base + ":" + bin;
      if (!buckets.has(key)) buckets.set(key, []);
      buckets.get(key).push(cluster);
    }
    cluster.observations.push(o);
    cluster.localRange = [
      Math.min(cluster.localRange[0], o.local.pos),
      Math.max(cluster.localRange[1], o.local.pos),
    ];
    cluster.partnerRange = [
      Math.min(cluster.partnerRange[0], o.partner.pos),
      Math.max(cluster.partnerRange[1], o.partner.pos),
    ];
    cluster.partner.start = Math.min(cluster.partner.start, o.partner.start);
    cluster.partner.end = Math.max(cluster.partner.end, o.partner.end);
  }
  if (invalidSA)
    warnings.push(
      `${invalidSA} malformed or oversized SA entries were ignored.`,
    );
  const byID = new Map(reads.map((r) => [r.id, r]));
  const events = clusters
    .map((c) => ({
      ...c,
      fragments: new Set(c.observations.map((o) => o.fragment)).size,
      eligibleFragments: new Set(
        c.observations.filter((o) => o.eligible).map((o) => o.fragment),
      ).size,
      localQualityFragments: new Set(
        c.observations.filter((o) => o.localOK).map((o) => o.fragment),
      ).size,
      unknownPartnerQuality: new Set(
        c.observations
          .filter((o) => o.partnerMapq === undefined || o.partnerMapq === 255)
          .map((o) => o.fragment),
      ).size,
      duplicateFragments: new Set(
        c.observations.filter((o) => o.duplicate).map((o) => o.fragment),
      ).size,
      records: [...new Set(c.observations.map((o) => o.record))],
      notes: [...new Set(c.observations.map((o) => o.note))],
      sources: [...new Set(c.observations.map((o) => o.source))],
      observations: undefined,
    }))
    .filter((e) => {
      // Structural links are clues rather than resolved alleles. Apply the same
      // threshold to their local-fragment evidence fraction, without querying a
      // partner locus or claiming a cancer-cell fraction.
      const covered = new Set(
        e.records
          .map((id) => byID.get(id))
          .filter(Boolean)
          .filter(
            (r) =>
              !(r.flags & (4 | 256 | 512)) &&
              (opts.duplicates || !(r.flags & 1024)) &&
              r.mq !== 255 &&
              r.mq >= opts.mapq,
          )
          .map(identity),
      );
      for (const r of reads) {
        if (
          (r.chr || region.chr) !== region.chr ||
          r.start > e.anchorLocal ||
          r.end <= e.anchorLocal ||
          r.flags & (4 | 256 | 512) ||
          (!opts.duplicates && r.flags & 1024) ||
          r.mq === 255 ||
          r.mq < opts.mapq ||
          (opts.sample && r.sample !== opts.sample)
        )
          continue;
        covered.add(identity(r));
      }
      e.localCallableFragments = covered.size;
      e.localEvidenceFraction = covered.size
        ? e.localQualityFragments / covered.size
        : null;
      return (
        e.localEvidenceFraction === null ||
        e.localEvidenceFraction >= (opts.minVaf ?? 0.02)
      );
    });
  const fingerprint = new Map();
  for (const e of events)
    for (const id of e.records) {
      if (!fingerprint.has(id)) fingerprint.set(id, []);
      fingerprint.get(id).push(e.id);
    }
  for (const r of reads) r.structuralLinks = fingerprint.get(r.id) || [];
  return {
    events,
    warnings,
    clusterTolerance: 100,
    distantMateThreshold: 10000,
  };
}
