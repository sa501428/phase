import { resolveVariantLines, requestedPairs } from "../src/variant-input.js";
import { phase } from "../src/analysis.js";

export function queryVariants(text, result, catalog = [], selected) {
  if (result.reference) {
    const entries = resolveVariantLines(text, {
      catalog,
      selected,
      region: result.region,
      reference: result.reference,
    });
    const called = phase(
      result.reads,
      result.reference,
      result.region.start,
      result.settings,
      entries.filter((e) => e.variant).map((e) => e.variant),
    );
    return {
      entries,
      pairs: requestedPairs(entries, called, {
        diploid: false,
        support: result.settings.support,
      }),
    };
  }
  const lines = text
    .split(/\r?\n/)
    .map((s) => s.trim())
    .filter(Boolean);
  if (lines.length > 32 || text.length > 32000)
    throw Error("Enter at most 32 variants / 32,000 characters");
  const entries = lines.map((input, i) => {
    let m = input
      .replace(/\s+/g, "")
      .match(/^([^:]+):(?:g\.)?(\d+)([ACGT])>([ACGT])$/i);
    if (!m) {
      const f = input.split(/\s+/);
      if (
        f.length === 4 &&
        /^\d+$/.test(f[1]) &&
        /^[ACGT]$/i.test(f[2]) &&
        /^[ACGT]$/i.test(f[3])
      )
        m = [null, ...f];
    }
    if (!m)
      return {
        input,
        line: i + 1,
        error:
          "Without a reference use chr:positionREF>ALT for SNVs. Add a reference for HGVS and indels.",
      };
    const pos = +m[2] - 1,
      ref = m[3].toUpperCase(),
      alt = m[4].toUpperCase();
    if (
      m[1].replace(/^chr/i, "") !== result.region.chr.replace(/^chr/i, "") ||
      pos < result.region.start ||
      pos >= result.region.end ||
      !Number.isSafeInteger(pos) ||
      alt === ref
    )
      return {
        input,
        line: i + 1,
        error: "Check contig, interval and distinct REF/ALT alleles",
      };
    return {
      input,
      line: i + 1,
      pos,
      ref,
      alt,
      variant: { pos, ref, alt },
      unverified: true,
    };
  });
  const pairs = [],
    resolved = entries.filter((e) => e.variant);
  for (let a = 0; a < resolved.length; a++)
    for (let b = a + 1; b < resolved.length; b++) {
      const left = resolved[a],
        right = resolved[b];
      let both = 0,
        leftOnly = 0,
        rightOnly = 0,
        neither = 0,
        unknown = 0,
        other = 0;
      for (const f of result.fragments) {
        if (f.conflict) {
          other++;
          continue;
        }
        const calls = new Map(f.calls),
          x = calls.get(left.pos),
          y = calls.get(right.pos);
        if (!x || !y) {
          unknown++;
          continue;
        }
        if (x === left.alt && y === right.alt) both++;
        else if (x === left.alt && y === right.ref) leftOnly++;
        else if (x === left.ref && y === right.alt) rightOnly++;
        else if (x === left.ref && y === right.ref) neither++;
        else other++;
      }
      pairs.push({
        left: left.input,
        right: right.input,
        both,
        leftOnly,
        rightOnly,
        neither,
        unknown,
        other,
        status:
          both >= result.settings.support
            ? leftOnly + rightOnly >= result.settings.support
              ? "Mixed patterns; cis also observed"
              : "Cis observed"
            : leftOnly >= result.settings.support &&
                rightOnly >= result.settings.support
              ? "Separate patterns; trans is not established"
              : "Unphased / insufficient shared evidence",
      });
    }
  return { entries, pairs };
}
