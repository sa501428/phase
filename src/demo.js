export function synthetic() {
  let seed = 7;
  const reference = Array.from({ length: 800 }, () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return "ACGT"[seed >>> 30];
  }).join("");
  const reads = [];
  const origin = 1000;
  for (let group = 0; group < 3; group++)
    for (let k = 0; k < 12; k++) {
      const start = 1050 + k * 2;
      let seq = reference.slice(start - origin, start - origin + 230),
        cigar = "230M";
      if (group > 0) {
        let a = seq.split("");
        a[1120 - start] = a[1120 - start] === "A" ? "T" : "A";
        a[1160 - start] = a[1160 - start] === "C" ? "G" : "C";
        seq = a.join("");
      }
      if (group === 2) {
        const at = 1200 - start,
          ins = reference.slice(190, 200);
        seq = seq.slice(0, at) + ins + seq.slice(at);
        cigar = `${at}M10I${230 - at}M`;
      }
      reads.push({
        id: reads.length,
        name: `synthetic-${group}-${k}`,
        start,
        end: start + 230,
        seq,
        cigar,
        qual: Array(seq.length).fill(35),
        flags: k % 2 ? 16 : 0,
        mq: 60,
        rg: "synthetic",
        sample: "synthetic",
        mate: -1,
        mateRef: -1,
      });
    }
  for (let k = 0; k < 8; k++) {
    const start = 1480 + k;
    const at = 1550 - start;
    const seq =
      reference.slice(start - origin, 550) +
      reference.slice(558, start - origin + 160);
    reads.push({
      id: reads.length,
      name: `deletion-${k}`,
      start,
      end: start + 160,
      seq,
      cigar: `${at}M8D${160 - at - 8}M`,
      qual: Array(seq.length).fill(35),
      flags: 0,
      mq: 60,
      rg: "synthetic",
      sample: "synthetic",
      mate: -1,
    });
  }
  return {
    reference,
    reads,
    region: { chr: "synthetic", start: origin, end: origin + 800 },
  };
}
