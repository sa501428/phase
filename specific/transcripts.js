export const reverseComplement = (s) =>
  [...s]
    .reverse()
    .map((b) => ({ A: "T", T: "A", C: "G", G: "C" })[b] || b)
    .join("");
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
