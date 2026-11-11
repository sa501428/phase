import {
  codingPosition,
  findTranscript,
  reverseComplement,
} from "./transcripts.js";

// All internal positions are zero-based; insertion positions are boundaries.
export const contigKey = (s) =>
  String(s).replace(/^chr/i, "").replace(/^MT$/i, "M").toUpperCase();
const dna = /^[ACGT]+$/i;
const aliases = {
  chr: ["chr", "chrom", "chromosome", "chromosome_name"],
  pos: ["pos", "position", "start", "start_position", "startposition"],
  ref: ["ref", "reference", "reference_allele", "referenceallele"],
  alt: [
    "alt",
    "alternate",
    "alternate_allele",
    "tumor_seq_allele2",
    "tumorseqallele2",
  ],
};
function columns(line) {
  if (line.includes("\t")) return line.split("\t").map((s) => s.trim());
  if (line.includes(",")) {
    const cells = [];
    let value = "",
      quoted = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (c === '"') {
        if (quoted && line[i + 1] === '"') {
          value += '"';
          i++;
        } else quoted = !quoted;
      } else if (c === "," && !quoted) {
        cells.push(value.trim());
        value = "";
      } else value += c;
    }
    if (quoted) throw Error("Unclosed CSV quote");
    cells.push(value.trim());
    return cells;
  }
  return line.trim().split(/\s+/);
}
function headerMap(fields) {
  const names = fields.map((s) => s.replace(/^#/, "").toLowerCase());
  const map = Object.fromEntries(
    Object.entries(aliases).map(([key, options]) => [
      key,
      names.findIndex((s) => options.includes(s)),
    ]),
  );
  return Object.values(map).every((i) => i >= 0) ? map : null;
}
function genomicContig(token, assembly) {
  const m = token.match(/^NC_0*(\d+)\.(\d+)$/i);
  if (!m) return token;
  const n = +m[1];
  const versions =
    assembly === "hg19"
      ? [
          10, 11, 11, 11, 9, 11, 13, 10, 11, 10, 9, 11, 10, 8, 9, 9, 10, 9, 9,
          10, 8, 10, 10, 9,
        ]
      : [
          11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10,
          10, 11, 9, 11, 11, 10,
        ];
  if (n === 12920 && +m[2] === 1) return "chrM";
  if (n < 1 || n > 24 || +m[2] !== versions[n - 1])
    throw Error(
      "Genomic accession does not match the selected assembly; use an explicit BAM contig",
    );
  return "chr" + (n === 23 ? "X" : n === 24 ? "Y" : n);
}
export function minimalVariant(v) {
  let { pos, ref, alt, refLength } = v;
  if (ref !== null && alt !== null) {
    while (ref.length && alt.length && ref.at(-1) === alt.at(-1)) {
      ref = ref.slice(0, -1);
      alt = alt.slice(0, -1);
    }
    while (ref.length && alt.length && ref[0] === alt[0]) {
      ref = ref.slice(1);
      alt = alt.slice(1);
      pos++;
    }
    refLength = ref.length;
  }
  if (
    !Number.isSafeInteger(pos) ||
    pos < 0 ||
    !Number.isSafeInteger(refLength) ||
    refLength < 0
  )
    throw Error("Coordinates must be positive safe integers");
  if (ref === alt && ref !== null)
    throw Error("Reference and alternate alleles are identical");
  if (refLength + (alt?.length || v.copyLength || 0) > 10000)
    throw Error(
      "Allele exceeds 10 kb; structural or uncertain events are not supported",
    );
  const result = { ...v, pos, ref, alt, refLength };
  result.key = JSON.stringify([
    contigKey(v.chr),
    pos,
    ref,
    alt,
    refLength,
    v.copyStart ?? null,
    v.copyLength ?? null,
  ]);
  result.label =
    ref === null
      ? `${v.chr}:g.${pos + 1}${refLength > 1 ? "_" + (pos + refLength) : ""}del${alt ? "ins" + alt : ""}`
      : alt === null
        ? `${v.chr}:g.${v.copyStart + 1}_${v.copyStart + v.copyLength}dup`
        : ref.length === 0
          ? `${v.chr}:g.${pos}_${pos + 1}ins${alt}`
          : alt.length === 0
            ? `${v.chr}:g.${pos + 1}${ref.length > 1 ? "_" + (pos + ref.length) : ""}del${ref}`
            : `${v.chr}:${pos + 1}${ref}>${alt}`;
  return result;
}
function fromAlleles(chr, position, ref, alternate) {
  const pos = Number(position) - 1;
  ref = ref.toUpperCase();
  alternate = alternate.toUpperCase();
  if (
    (!dna.test(ref) && ref !== "-") ||
    (!dna.test(alternate) && alternate !== "-")
  )
    throw Error(
      "Use explicit A/C/G/T alleles; symbolic, breakend, missing and star alleles are not supported",
    );
  return minimalVariant({
    chr,
    pos,
    ref: ref === "-" ? "" : ref,
    alt: alternate === "-" ? "" : alternate,
    refLength: ref === "-" ? 0 : ref.length,
  });
}
function fromHGVS(identifier, type, expression, options) {
  type = type.toLowerCase();
  const context =
    type === "c" ? findTranscript(identifier, options.catalog || []) : null;
  const chr = context
    ? context.chr
    : genomicContig(identifier || "", options.assembly);
  if (!chr) throw Error("Include a contig for g. or an accession for c.");
  const position = "(?:\\*\\d+|-?[1-9]\\d*)(?:[+-]\\d+)?";
  const e = expression.match(
    new RegExp(
      `^(${position})(?:_(${position}))?(?:([ACGT]+)>([ACGT]+)|delins([ACGT]+)|(del|dup)([ACGT]*)|ins([ACGT]+))$`,
      "i",
    ),
  );
  if (!e)
    throw Error(
      "Unsupported HGVS; use a DNA substitution, del, ins, dup or delins with exact coordinates",
    );
  const coord = (s) =>
    context
      ? codingPosition(s, context)
      : /^[1-9]\d*$/.test(s)
        ? Number(s) - 1
        : NaN;
  const a = coord(e[1]),
    b = coord(e[2] || e[1]),
    minus = context?.strand === "-";
  if (!Number.isSafeInteger(a) || !Number.isSafeInteger(b) || a < 0 || b < 0)
    throw Error("Invalid genomic position");
  if (minus ? a < b : a > b)
    throw Error("Range must follow the genomic/transcript direction");
  if (context && e[2]) {
    const exons = context.transcript[3];
    const first = exons.findIndex(([s, end]) => a >= s && a < end),
      last = exons.findIndex(([s, end]) => b >= s && b < end);
    if (first >= 0 && last >= 0 && first !== last)
      throw Error(
        "Range crosses a splice junction; enter the genomic DNA event",
      );
  }
  const start = Math.min(a, b),
    end = Math.max(a, b) + 1;
  const orient = (s) =>
    minus ? reverseComplement(s.toUpperCase()) : s.toUpperCase();
  let pos = start,
    ref = null,
    alt = "",
    refLength = end - start,
    copyStart,
    copyLength;
  if (e[3]) {
    ref = orient(e[3]);
    alt = orient(e[4]);
    refLength = ref.length;
    if (e[2] && refLength !== end - start)
      throw Error("Substitution range and REF length differ");
    if (context && !e[2] && refLength !== 1)
      throw Error("Use c. delins with a range for multi-base substitutions");
  } else if (e[8]) {
    if (!e[2] || Math.abs(a - b) !== 1)
      throw Error("Insertion requires adjacent flanking positions");
    pos = Math.max(a, b);
    ref = "";
    refLength = 0;
    alt = orient(e[8]);
  } else if (e[5]) alt = orient(e[5]);
  else {
    const specified = e[7] ? orient(e[7]) : null;
    if (specified && specified.length !== refLength)
      throw Error(
        "Deleted/duplicated sequence length does not match the range",
      );
    if (e[6].toLowerCase() === "dup") {
      pos = minus ? start : end;
      ref = "";
      refLength = 0;
      alt = specified;
      if (!specified) {
        copyStart = start;
        copyLength = end - start;
      }
    } else ref = specified;
  }
  return minimalVariant({
    chr,
    pos,
    ref,
    alt,
    refLength,
    copyStart,
    copyLength,
    transcript: context?.transcript[0],
  });
}
export function parseVariants(text, options = {}) {
  const entries = [],
    variants = [],
    seen = new Map();
  let header = null;
  const add = (variant, line, input) => {
    let id = seen.get(variant.key);
    const duplicate = id !== undefined;
    if (!duplicate) {
      id = variants.length;
      seen.set(variant.key, id);
      variants.push({ ...variant, id, sources: [] });
    }
    variants[id].sources.push(line);
    entries.push({ line, input, id, label: variants[id].label, duplicate });
  };
  for (const [index, original] of text.split(/\r?\n/).entries()) {
    const input = original.trim(),
      line = index + 1;
    if (!input || input.startsWith("##")) continue;
    try {
      const fields = columns(input),
        possibleHeader = headerMap(fields);
      if (possibleHeader) {
        header = possibleHeader;
        continue;
      }
      if (input.startsWith("#")) continue;
      if (header) {
        const chr = fields[header.chr],
          pos = fields[header.pos],
          ref = fields[header.ref],
          alt = fields[header.alt];
        if (!chr || !/^\d+$/.test(pos || "") || !ref || !alt)
          throw Error("Row is missing chromosome, position, REF or ALT");
        for (const a of alt.split(","))
          add(fromAlleles(chr, pos, ref, a), line, input);
        continue;
      }
      let m;
      const coding = input.match(/:c\.([^:\s);,]+)/i);
      const accessions = [
        ...new Set(input.match(/(?:NM|NR|XM|XR)_\d+(?:\.\d+)?/gi) || []),
      ];
      if (coding && accessions.length) {
        if (accessions.length > 1)
          throw Error(
            "Multiple transcript accessions on one line; use one transcript annotation per line",
          );
        const expression = coding[1].replace(
          /^([ACGT])(\d+)([ACGT])$/i,
          "$2$1>$3",
        );
        add(fromHGVS(accessions[0], "c", expression, options), line, input);
      } else if (
        (m = input.match(/(?:^|[\s(;])([^\s():;]+):([cg])\.([^\s);,]+)/i))
      ) {
        add(fromHGVS(m[1], m[2], m[3], options), line, input);
      } else if (
        (m = input.match(
          /(?:^|[\s(;])([^\s():;]+):(?:g\.)?(\d+)([ACGT]+)>([ACGT]+)(?=$|[\s);,])/i,
        ))
      ) {
        add(fromAlleles(m[1], m[2], m[3], m[4]), line, input);
      } else if (
        (m = input.match(
          /^([^:\s]+)[:\-](\d+)[:\-]([ACGT-]+)[:\-]([ACGT,\-]+)$/i,
        ))
      ) {
        for (const a of m[4].split(","))
          add(fromAlleles(m[1], m[2], m[3], a), line, input);
      } else if (fields.length >= 4 && /^\d+$/.test(fields[1])) {
        // Standard VCF has ID in column 3; a four-column table does not.
        const vcf =
          fields.length >= 5 && (dna.test(fields[3]) || fields[3] === "-");
        const ref = fields[vcf ? 3 : 2],
          alt = fields[vcf ? 4 : 3];
        for (const a of alt.split(","))
          add(fromAlleles(fields[0], fields[1], ref, a), line, input);
      } else
        throw Error(
          "Use CHR POS REF ALT, VCF/MAF/CSV with a header, chr:posREF>ALT, chr:pos:REF:ALT, or accession:c./chr:g. HGVS. Protein-only and rs IDs need genomic alleles.",
        );
    } catch (error) {
      entries.push({ line, input, error: error.message });
    }
  }
  return { entries, variants };
}

export function resolveContigs(variants, contigs) {
  return variants.map((v) => {
    const exact = contigs.find((c) => c.refName === v.chr);
    const matches = exact
      ? [exact]
      : contigs.filter((c) => contigKey(c.refName) === contigKey(v.chr));
    if (matches.length !== 1)
      throw Error(`${v.label}: BAM contig is missing or ambiguous`);
    const c = matches[0];
    const end = Math.max(
      v.pos + Math.max(v.refLength, 1),
      (v.copyStart || 0) + (v.copyLength || 0),
    );
    if (
      v.pos < 0 ||
      end > c.length ||
      (!v.refLength && (v.pos < 1 || v.pos >= c.length))
    )
      throw Error(
        `${v.label}: event or its insertion flanks fall outside the BAM contig`,
      );
    return { ...v, chr: c.refName };
  });
}
