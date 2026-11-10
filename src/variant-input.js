import { normalize } from "./analysis.js";
export const reverseComplement = (s) =>
  [...s]
    .reverse()
    .map((b) => ({ A: "T", T: "A", C: "G", G: "C" })[b] || b)
    .join("");
const chrKey = (s) => s.replace(/^chr/i, "").replace(/^MT$/, "M");
export function findTranscript(identifier, catalog, selected) {
  if (!identifier) {
    if (!selected) throw Error("Select a transcript for a bare c. entry");
    return selected;
  }
  const choices = [];
  for (const g of catalog)
    for (const t of g[6])
      if (
        t[0] === identifier ||
        (!identifier.includes(".") && t[0].split(".")[0] === identifier) ||
        g[0].toUpperCase() === identifier.toUpperCase()
      )
        choices.push({ chr: g[1], strand: g[2], transcript: t });
  if (
    selected &&
    choices.some(
      (c) =>
        c.chr === selected.chr && c.transcript[0] === selected.transcript[0],
    )
  )
    return selected;
  if (choices.length !== 1)
    throw Error(
      choices.length
        ? "Ambiguous transcript; specify accession.version or select one"
        : "Transcript not in the offline catalog for this assembly",
    );
  return choices[0];
}
export function codingPosition(token, context) {
  const t = context.transcript,
    exons = [...t[3]].sort((a, b) => a[0] - b[0]);
  if (context.strand === "-") exons.reverse();
  const direction = context.strand === "-" ? -1 : 1;
  if (t[6] !== undefined && (t[6] !== "cmpl" || t[7] !== "cmpl"))
    throw Error(
      "Transcript CDS is incomplete; use an explicit genomic variant",
    );
  if (t[4] === undefined || t[5] === undefined || t[4] >= t[5])
    throw Error("Transcript has no coding sequence");
  function offsetOf(genomic) {
    let n = 0;
    for (const [s, e] of exons) {
      if (genomic >= s && genomic < e)
        return n + (direction === 1 ? genomic - s : e - 1 - genomic);
      n += e - s;
    }
    throw Error("CDS boundary outside transcript exons");
  }
  const start = offsetOf(direction === 1 ? t[4] : t[5] - 1),
    end = offsetOf(direction === 1 ? t[5] - 1 : t[4]);
  const m = token.match(/^(\*\d+|-?[1-9]\d*)([+-]\d+)?$/);
  if (!m) throw Error("Invalid c. position (c.0 is not defined)");
  if (m[1] === "*0" || (m[2] && Number(m[2]) === 0))
    throw Error("HGVS zero coordinates/offsets are not defined");
  const raw = m[1],
    n = raw.startsWith("*")
      ? end + Number(raw.slice(1))
      : Number(raw) < 0
        ? start + Number(raw)
        : start + Number(raw) - 1;
  if ((!raw.startsWith("*") && Number(raw) > 0 && n > end) || n < 0)
    throw Error("Coding position outside CDS/UTR; use * for the 3-prime UTR");
  let left = n;
  for (let i = 0; i < exons.length; i++) {
    const [s, e] = exons[i],
      length = e - s;
    if (left < length) {
      let pos = direction === 1 ? s + left : e - 1 - left;
      if (m[2]) {
        const off = Number(m[2]),
          plus = off > 0;
        if ((plus && left !== length - 1) || (!plus && left !== 0))
          throw Error("Intronic offset must be anchored at an exon boundary");
        const next = exons[i + (plus ? 1 : -1)];
        if (!next) throw Error("Intronic offset has no adjacent exon");
        const target = pos + direction * off,
          min = Math.min(
            pos,
            plus
              ? direction === 1
                ? next[0]
                : next[1] - 1
              : direction === 1
                ? next[1] - 1
                : next[0],
          ),
          max = Math.max(
            pos,
            plus
              ? direction === 1
                ? next[0]
                : next[1] - 1
              : direction === 1
                ? next[1] - 1
                : next[0],
          );
        if (target <= min || target >= max)
          throw Error("Intronic offset reaches an exon or exceeds the intron");
        pos = target;
      }
      return pos;
    }
    left -= length;
  }
  throw Error("Coding position outside transcript");
}
export function resolveVariantLines(
  text,
  { catalog = [], selected, region, reference },
) {
  if (text.length > 32000)
    throw Error("Variant entry exceeds 32,000 characters");
  const lines = text
    .split(/\r?\n/)
    .map((s) => s.trim())
    .filter(Boolean);
  if (lines.length > 32) throw Error("Enter at most 32 variants at a time");
  return lines.map((input, i) => {
    try {
      let chr,
        expression,
        context = null,
        type = "g";
      const compact = input.replace(/\s+/g, "");
      let m = compact.match(/^(?:(.+):)?([cg])\.(.+)$/i);
      if (m) {
        type = m[2].toLowerCase();
        expression = m[3];
        if (type === "c") {
          context = findTranscript(m[1], catalog, selected);
          chr = context.chr;
        } else chr = m[1] || region.chr;
      } else if ((m = compact.match(/^([^:]+):(\d+)([ACGT]+)>([ACGT]+)$/i))) {
        chr = m[1];
        expression = m[2] + m[3] + ">" + m[4];
      } else {
        const fields = input.split(/[\s,:]+/);
        if (
          fields.length === 4 &&
          /^\d+$/.test(fields[1]) &&
          /^[ACGT]+$/i.test(fields[2] + fields[3])
        ) {
          chr = fields[0];
          expression = fields[1] + fields[2] + ">" + fields[3];
        } else
          throw Error(
            "Use chr:positionREF>ALT, chr:g.HGVS, transcript:c.HGVS, or CHR POS REF ALT",
          );
      }
      if (chrKey(chr) !== chrKey(region.chr))
        throw Error(
          "Variant is on another contig; analyze that locus separately",
        );
      const posToken = "(?:\\*\\d+|-?[1-9]\\d*)(?:[+-]\\d+)?";
      const re = new RegExp(
        "^(" +
          posToken +
          ")(?:_(" +
          posToken +
          "))?(?:([ACGT]+)>([ACGT]+)|(delins)([ACGT]+)|(del|dup)([ACGT]*)|(ins)([ACGT]+))$",
        "i",
      );
      m = expression.match(re);
      if (!m)
        throw Error(
          "Unsupported or uncertain HGVS; use a DNA substitution, del, ins, dup, or delins",
        );
      const coord = (s) =>
        type === "c"
          ? codingPosition(s, context)
          : /^\d+$/.test(s)
            ? Number(s) - 1
            : (() => {
                throw Error("Genomic positions must be positive integers");
              })();
      const a = coord(m[1]),
        b = coord(m[2] || m[1]),
        minus = context?.strand === "-",
        start = Math.min(a, b),
        end = Math.max(a, b) + 1;
      if ((!minus && a > b) || (minus && a < b))
        throw Error("Range must follow reference/transcript direction");
      const sequence = (s, e) => {
        if (s < region.start || e > region.end)
          throw Error(
            "Outside the analyzed interval; enlarge the region and rerun",
          );
        return reference.slice(s - region.start, e - region.start);
      };
      if (type === "c" && m[2]) {
        const exons = context.transcript[3],
          first = exons.findIndex(([s, e]) => a >= s && a < e),
          last = exons.findIndex(([s, e]) => b >= s && b < e);
        if (first >= 0 && last >= 0 && first !== last)
          throw Error(
            "Range crosses a splice junction; enter the genomic DNA event explicitly",
          );
      }
      let pos = start,
        ref = "",
        alt = "";
      if (m[3]) {
        ref = m[3].toUpperCase();
        alt = m[4].toUpperCase();
        if (m[2] || ref.length !== 1 || alt.length !== 1) {
          if (!m[2]) {
            if (type === "c")
              throw Error("Use delins for multi-base coding substitutions");
            ref =
              sequence(start, start + ref.length) === ref
                ? ref
                : (() => {
                    throw Error("Reference allele mismatch");
                  })();
          } else if (end - start !== ref.length)
            throw Error("Substitution range/allele length mismatch");
        }
        if (minus) {
          ref = reverseComplement(ref);
          alt = reverseComplement(alt);
        }
      } else if (m[9]) {
        if (!m[2] || Math.abs(a - b) !== 1)
          throw Error(
            "Insertion requires two adjacent genomic flanking positions",
          );
        pos = Math.max(a, b);
        sequence(pos - 1, pos + 1);
        alt = m[10].toUpperCase();
        if (minus) alt = reverseComplement(alt);
      } else {
        ref = sequence(start, end);
        const op = (m[5] || m[7]).toLowerCase(),
          specified = m[8]?.toUpperCase();
        if (
          specified &&
          (minus ? reverseComplement(specified) : specified) !== ref
        )
          throw Error("Deleted/duplicated sequence does not match reference");
        if (op === "delins") {
          alt = m[6].toUpperCase();
          if (minus) alt = reverseComplement(alt);
        } else if (op === "dup") {
          alt = ref;
          ref = "";
          pos = minus ? start : end;
        }
      }
      if (ref && sequence(pos, pos + ref.length) !== ref)
        throw Error("Reference allele mismatch: check assembly and transcript");
      if (ref === alt)
        throw Error("Reference and alternate alleles are identical");
      if (ref.length + alt.length > 5000)
        throw Error("Requested allele exceeds 5 kb limit");
      const variant = {
        ...normalize(pos, ref, alt, reference, region.start),
        supplied: true,
      };
      return {
        line: i + 1,
        input,
        chr: region.chr,
        pos,
        ref,
        alt,
        transcript: context?.transcript[0],
        variant,
      };
    } catch (error) {
      return { line: i + 1, input, error: error.message };
    }
  });
}
export function requestedPairs(
  entries,
  result,
  { diploid = false, support = 3 } = {},
) {
  const resolved = entries.filter((e) => e.variant),
    pairs = [];
  for (let a = 0; a < resolved.length; a++)
    for (let b = a + 1; b < resolved.length; b++) {
      const left = resolved[a],
        right = resolved[b],
        i = result.variants.findIndex((v) => v.id === left.variant.id),
        j = result.variants.findIndex((v) => v.id === right.variant.id);
      let both = 0,
        leftOnly = 0,
        rightOnly = 0,
        neither = 0,
        other = 0,
        unknown = 0;
      for (const v of result.vectors) {
        const x = v[i],
          y = v[j];
        if (x === undefined || y === undefined || x === "?" || y === "?") {
          unknown++;
          continue;
        }
        const xa = x === (left.variant.alt || "-"),
          ya = y === (right.variant.alt || "-"),
          xr = x === (left.variant.ref || "-"),
          yr = y === (right.variant.ref || "-");
        if (xa && ya) both++;
        else if (xa && yr) leftOnly++;
        else if (xr && ya) rightOnly++;
        else if (xr && yr) neither++;
        else other++;
      }
      let status = "Unphased";
      if (left.variant.id === right.variant.id) status = "Same variant";
      else if (both && (leftOnly || rightOnly)) status = "Mixed patterns";
      else if (both >= support) status = "Cis supported";
      else if (both) status = "Possible cis";
      else if (leftOnly >= support && rightOnly >= support)
        status = diploid
          ? "Trans supported (diploid assumption)"
          : "Separate patterns; trans unresolved";
      pairs.push({
        a: left.line,
        b: right.line,
        status,
        both,
        leftOnly,
        rightOnly,
        neither,
        other,
        unknown,
        informative: both + leftOnly + rightOnly + neither + other,
      });
    }
  return pairs;
}
