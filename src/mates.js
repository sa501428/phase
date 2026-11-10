// Fetch named, mapped mates on this contig using the BAM index, not a file scan.
export async function completePairs(bam, records, contig, refID) {
  const result = [...records],
    seen = new Set(records.map((r) => r.fileOffset)),
    present = new Set(
      records.map((r) =>
        JSON.stringify([
          r.name,
          String(r.getTag("RG") || ""),
          r.start,
          r.flags & 192,
        ]),
      ),
    ),
    queries = new Map();
  for (const r of records) {
    if (
      !(r.flags & 1) ||
      r.flags & 8 ||
      r.next_refid !== refID ||
      r.next_pos < 0 ||
      r.flags & (256 | 2048)
    )
      continue;
    const opposite = r.flags & 64 ? 128 : r.flags & 128 ? 64 : 0;
    if (!opposite) continue;
    const key = JSON.stringify([
      r.name,
      String(r.getTag("RG") || ""),
      r.next_pos,
      opposite,
    ]);
    if (present.has(key)) continue;
    const bin = Math.floor(r.next_pos / 1000) * 1000;
    if (!queries.has(bin)) queries.set(bin, new Set());
    queries.get(bin).add(key);
  }
  if (queries.size > 512)
    throw Error(
      "Mate lookup exceeds 512 indexed windows; narrow the interval. No reads were downsampled.",
    );
  for (const [start, wanted] of queries) {
    for (const r of await bam.getRecordsForRange(contig, start, start + 1000)) {
      const key = JSON.stringify([
        r.name,
        String(r.getTag("RG") || ""),
        r.start,
        r.flags & 192,
      ]);
      if (wanted.has(key) && !seen.has(r.fileOffset)) {
        seen.add(r.fileOffset);
        result.push(r);
      }
    }
    if (result.length > 100000)
      throw Error(
        "Mate retrieval exceeds 100,000 records; narrow the interval",
      );
  }
  return result;
}
