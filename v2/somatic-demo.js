// Ground truth describes sequence lineages used to generate the mixture. It
// is never supplied to the assembler or used as a requested clone count.
export function somaticDemo({ partials = false } = {}) {
  let seed = 19;
  const reference = Array.from({ length: 800 }, () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return "ACGT"[seed >>> 30];
  }).join("");
  const mutations = [
    [1050],
    [1050, 1130],
    [1050, 1170],
    [1050, 1210],
    [1050, 1130, 1410],
    [1090],
    [1090, 1250],
    [1090, 1290],
    [1090, 1330],
    [1090, 1290, 1450],
  ];
  const abundances = [36, 24, 18, 12, 6, 30, 21, 15, 9, 3],
    reads = [];
  const alt = (p) => (reference[p - 1000] === "A" ? "T" : "A");
  function record(start, length, mutations, name, flags, mate) {
    let seq = reference.slice(start - 1000, start - 1000 + length).split("");
    for (const pos of mutations)
      if (pos >= start && pos < start + length) seq[pos - start] = alt(pos);
    seq = seq.join("");
    return {
      id: reads.length,
      name,
      start,
      end: start + length,
      seq,
      cigar: length + "M",
      qual: Array(length).fill(35),
      flags,
      mq: 60,
      chr: "synthetic",
      mate,
      mateChr: "synthetic",
      mateMq: 60,
      rg: "synthetic",
      sample: "synthetic",
    };
  }
  for (let h = 0; h < mutations.length; h++)
    for (let n = 0; n < abundances[h]; n++) {
      const start = 1020 + (n % 5) * 2,
        mate = start + 250,
        name = `lineage-${h + 1}-fragment-${n + 1}`;
      reads.push(record(start, 300, mutations[h], name, 99, mate));
      reads.push(record(mate, 300, mutations[h], name, 147, start));
    }
  if (partials)
    for (let n = 0; n < 5; n++)
      reads.push(record(1015, 100, [1050], `shared-founder-A-${n}`, 0, -1));
  return {
    reads,
    reference,
    region: { chr: "synthetic", start: 1000, end: 1800 },
    truth: { mutations, abundances, alt },
  };
}
