export function cigarOps(cigar) {
  const out = [...cigar.matchAll(/(\d+)([MIDNSHP=X])/g)].map((x) => [
    +x[1],
    x[2],
  ]);
  if (
    out.map((x) => x.join("")).join("") !== cigar ||
    out.some((x) => x[0] < 1 || x[0] > 1e7)
  )
    throw Error("Invalid CIGAR");
  return out;
}
