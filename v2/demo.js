export function demo() {
  let seed = 7;
  const reference = Array.from({ length: 1000 }, () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return "ACGT"[seed >>> 30];
  }).join("");
  const origin = 1000,
    reads = [];
  function record(start, end, edits, extra) {
    let seq = "",
      cigar = [],
      p = start,
      run = 0;
    const flush = () => {
      if (run) cigar.push(run + "M");
      run = 0;
    };
    while (p < end) {
      const e = edits.find((e) => e.pos === p);
      if (e?.ins) {
        flush();
        seq += e.ins;
        cigar.push(e.ins.length + "I");
      }
      if (e?.del) {
        flush();
        cigar.push(e.del + "D");
        p += e.del;
        continue;
      }
      seq += e?.alt || reference[p - origin];
      p++;
      run++;
    }
    flush();
    return {
      id: reads.length,
      start,
      end,
      seq,
      cigar: cigar.join(""),
      qual: Array(seq.length).fill(35),
      mq: 60,
      rg: "demo",
      sample: "demo",
      chr: "synthetic",
      ...extra,
    };
  }
  for (let h = 0; h < 3; h++)
    for (let n = 0; n < 12; n++) {
      const edits = [];
      if (h > 0)
        edits.push(
          { pos: 1120, alt: "A" === reference[120] ? "T" : "A" },
          { pos: 1160, alt: "C" === reference[160] ? "G" : "C" },
        );
      if (h === 2) edits.push({ pos: 1200, ins: reference.slice(190, 200) });
      const start = 1040 + n * 2,
        mate = start + 70,
        name = `haplotype-${h + 1}-pair-${n + 1}`;
      reads.push(
        record(start, start + 220, edits, {
          name,
          flags: 99,
          mate,
          mateChr: "synthetic",
          mateMq: 60,
        }),
      );
      reads.push(
        record(mate, mate + 220, edits, {
          name,
          flags: 147,
          mate: start,
          mateChr: "synthetic",
          mateMq: 60,
        }),
      );
    }
  for (let h = 0; h < 2; h++)
    for (let n = 0; n < 8; n++) {
      const start = 1530 + n * 2,
        mate = start + 55,
        name = `deletion-${h}-pair-${n}`;
      const edits = h ? [{ pos: 1640, del: 8 }] : [];
      reads.push(
        record(start, start + 180, edits, {
          name,
          flags: 99,
          mate,
          mateChr: "synthetic",
          mateMq: 60,
        }),
      );
      reads.push(
        record(mate, mate + 180, edits, {
          name,
          flags: 147,
          mate: start,
          mateChr: "synthetic",
          mateMq: 60,
        }),
      );
    }
  return {
    reads,
    reference,
    region: { chr: "synthetic", start: origin, end: origin + 1000 },
  };
}
