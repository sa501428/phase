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
export function normalize(pos, ref, alt, sequence, origin) {
  while (ref.length && alt.length && ref.at(-1) === alt.at(-1)) {
    ref = ref.slice(0, -1);
    alt = alt.slice(0, -1);
  }
  while (ref.length && alt.length && ref[0] === alt[0]) {
    pos++;
    ref = ref.slice(1);
    alt = alt.slice(1);
  }
  if (!ref.length || !alt.length) {
    let allele = ref || alt;
    while (pos > origin && sequence[pos - origin - 1] === allele.at(-1)) {
      pos--;
      allele = allele.at(-1) + allele.slice(0, -1);
    }
    if (ref) ref = allele;
    else alt = allele;
  }
  return { pos, ref, alt, id: `${pos}:${ref}>${alt}` };
}
export function project(read, reference, origin, minBQ) {
  let p = read.start,
    q = 0;
  const bases = new Map(),
    events = [],
    clips = [],
    segments = [];
  for (const [n, op] of cigarOps(read.cigar)) {
    if ("M=X".includes(op)) {
      segments.push([p, p + n]);
      for (let j = 0; j < n; j++) {
        const b = read.seq[q + j],
          qual = read.qual?.[q + j];
        if (
          b &&
          "ACGT".includes(b) &&
          qual !== undefined &&
          qual !== 255 &&
          qual >= minBQ
        )
          bases.set(p + j, { allele: b, q: qual });
      }
      p += n;
      q += n;
    } else if (op === "I" || op === "D") {
      const ref = op === "D" ? reference.slice(p - origin, p - origin + n) : "",
        alt = op === "I" ? read.seq.slice(q, q + n) : "";
      const qual =
        op === "I"
          ? Math.min(...read.qual.slice(q, q + n))
          : Math.min(read.qual[q - 1] ?? 0, read.qual[q] ?? 0);
      if (
        p >= origin &&
        p + (op === "D" ? n : 0) <= origin + reference.length &&
        qual >= minBQ &&
        qual !== 255 &&
        /^[ACGT]*$/.test(ref + alt)
      )
        events.push({
          ...normalize(p, ref, alt, reference, origin),
          rawPos: p,
          q: qual,
        });
      if (op === "I") q += n;
      else p += n;
    } else if (op === "S") {
      clips.push({
        pos: p,
        side: q === 0 ? "left" : "right",
        seq: read.seq.slice(q, q + n),
        quality:
          read.qual.length >= q + n
            ? Math.min(...read.qual.slice(q, q + n))
            : 0,
      });
      q += n;
    } else if (op === "N") p += n;
  }
  if (q !== read.seq.length)
    throw Error("Sequence length does not match CIGAR");
  return { bases, events, clips, segments };
}
export function clipEvent(c, reference, origin) {
  const anchor = c.pos - origin;
  if (c.seq.length < 16 || c.seq.length > 1000) return null;
  // Require unique exact flanking sequence within a bounded local context.
  const lo = Math.max(0, anchor - 5000),
    hi = Math.min(reference.length, anchor + 5000),
    window = reference.slice(lo, hi);
  if (c.side === "right") {
    const prefix = c.seq.slice(0, 16),
      hit = window.indexOf(prefix);
    if (hit >= 0 && window.indexOf(prefix, hit + 1) < 0) {
      const target = lo + hit;
      if (
        target > anchor &&
        target - anchor <= 5000 &&
        reference.slice(target, target + c.seq.length) === c.seq
      )
        return {
          ...normalize(
            c.pos,
            reference.slice(anchor, target),
            "",
            reference,
            origin,
          ),
          q: 30,
          clipDerived: true,
        };
    }
    for (let n = c.seq.length - 12; n >= 1; n--) {
      const tail = c.seq.slice(n);
      if (
        reference.slice(anchor, anchor + tail.length) === tail &&
        tail.length >= 12 &&
        window.indexOf(tail) === window.lastIndexOf(tail)
      )
        return {
          ...normalize(c.pos, "", c.seq.slice(0, n), reference, origin),
          q: 30,
          clipDerived: true,
        };
    }
  } else {
    for (let n = 1; n <= c.seq.length - 12; n++) {
      const head = c.seq.slice(0, c.seq.length - n);
      if (
        reference.slice(anchor - head.length, anchor) === head &&
        head.length >= 12 &&
        window.indexOf(head) === window.lastIndexOf(head)
      )
        return {
          ...normalize(c.pos, "", c.seq.slice(-n), reference, origin),
          q: 30,
          clipDerived: true,
        };
    }
  }
  return null;
}
export function fragments(reads, reference, origin, opts) {
  const map = new Map();
  let excluded = 0;
  for (const r of reads) {
    r.excluded =
      !!(r.flags & (4 | 256 | 512)) ||
      (!opts.duplicates && !!(r.flags & 1024)) ||
      r.mq === 255 ||
      r.mq < opts.mapq ||
      (opts.sample && r.sample !== opts.sample);
    if (r.excluded) {
      excluded++;
      continue;
    }
    r.projection = project(r, reference, origin, opts.baseq);
    for (const clip of r.projection.clips) {
      if (clip.quality < opts.baseq || clip.quality === 255) continue;
      const event = clipEvent(clip, reference, origin);
      if (event) r.projection.events.push(event);
    }
    const key = JSON.stringify([
      r.sample,
      r.rg,
      r.name,
      r.flags & 1 || r.flags & 2048 ? "template" : r.id,
    ]);
    let f = map.get(key);
    if (!f) {
      f = {
        id: map.size,
        reads: [],
        bases: new Map(),
        events: new Map(),
        conflicts: new Set(),
      };
      map.set(key, f);
    }
    f.reads.push(r);
    for (const [p, b] of r.projection.bases) {
      const prev = f.bases.get(p);
      if (prev && prev.allele !== b.allele) {
        f.bases.delete(p);
        f.conflicts.add(p);
      } else if (!f.conflicts.has(p)) f.bases.set(p, b);
    }
    for (const e of r.projection.events) f.events.set(e.id, e);
  }
  return { items: [...map.values()], excluded };
}
export function parseVCF(text, contig, reference, origin) {
  const out = [];
  for (const line of text.split(/\r?\n/)) {
    if (!line || line[0] === "#") continue;
    const [chr, p, , ref, alts] = line.split("\t");
    if (chr.replace(/^chr/, "") !== contig.replace(/^chr/, "")) continue;
    const pos = +p - 1;
    if (pos < origin || pos >= origin + reference.length) continue;
    if (
      reference.slice(pos - origin, pos - origin + ref.length) !==
      ref.toUpperCase()
    )
      throw Error("VCF reference allele mismatch");
    for (const alt of alts.split(",")) {
      if (!/^[ACGT]+$/i.test(ref + alt)) continue;
      out.push({
        ...normalize(
          pos,
          ref.toUpperCase(),
          alt.toUpperCase(),
          reference,
          origin,
        ),
        supplied: true,
      });
    }
  }
  return out;
}
function allele(f, v) {
  if (!v.ref.length || !v.alt.length) {
    if (f.events.has(v.id)) {
      const end = v.pos + Math.max(1, v.ref.length);
      const contradiction = f.reads.some(
        (r) =>
          !r.projection.events.some((e) => e.id === v.id) &&
          r.projection.segments.some(([s, e]) => s < v.pos && e > end) &&
          r.projection.bases.has(v.pos - 1) &&
          r.projection.bases.has(end),
      );
      return contradiction ? "?" : v.alt || "-";
    }
    // absence of an indel is evidence only when a single alignment spans both anchors
    const end = v.pos + Math.max(1, v.ref.length);
    if (
      [...f.events.values()].some(
        (e) =>
          Math.abs(e.pos - v.pos) < Math.max(e.ref.length, v.ref.length, 1),
      )
    )
      return "?";
    return f.reads.some(
      (r) =>
        r.projection.segments.some(([s, e]) => s < v.pos && e > end) &&
        r.projection.bases.has(v.pos - 1) &&
        r.projection.bases.has(end),
    )
      ? v.ref || "-"
      : "?";
  }
  let s = "";
  for (let p = v.pos; p < v.pos + v.ref.length; p++) {
    const b = f.bases.get(p);
    if (!b) return "?";
    s += b.allele;
  }
  return s;
}
export function phase(reads, reference, origin, opts, supplied = []) {
  const { items: fs, excluded } = fragments(reads, reference, origin, opts),
    candidate = new Map(supplied.map((v) => [v.id, { ...v, support: 0 }]));
  const counts = new Map();
  for (const f of fs) {
    const seen = new Map(f.events);
    for (const [p, b] of f.bases) {
      const ref = reference[p - origin];
      if (ref && ref !== b.allele && "ACGT".includes(ref)) {
        const v = normalize(p, ref, b.allele, reference, origin);
        seen.set(v.id, v);
      }
    }
    for (const v of seen.values()) {
      const entry = counts.get(v.id) || { ...v, support: 0 };
      entry.support++;
      counts.set(v.id, entry);
    }
  }
  for (const v of counts.values())
    if (v.support >= opts.support || candidate.has(v.id))
      candidate.set(v.id, { ...v, supplied: candidate.has(v.id) });
  let variants = [...candidate.values()]
    .filter((v) => v.pos >= origin && v.pos < origin + reference.length)
    .sort((a, b) => a.pos - b.pos || a.id.localeCompare(b.id));
  if (variants.length > 500)
    throw Error(
      "More than 500 candidate alleles; narrow the interval or increase support",
    );
  variants = variants.map((v, i) => ({
    ...v,
    index: i,
    support: counts.get(v.id)?.support || 0,
    overlap: variants.some(
      (w) =>
        w.id !== v.id &&
        w.pos <= v.pos + Math.max(1, v.ref.length) - 1 &&
        v.pos <= w.pos + Math.max(1, w.ref.length) - 1,
    ),
  }));
  const vectors = fs.map((f) => variants.map((v) => allele(f, v))),
    edges = new Map(),
    parent = variants.map((_, i) => i);
  const root = (i) => (parent[i] === i ? i : (parent[i] = root(parent[i])));
  for (const a of vectors) {
    const idx = a.map((x, i) => (x !== "?" ? i : -1)).filter((i) => i >= 0);
    for (let j = 1; j < idx.length; j++) {
      const x = idx[j - 1],
        y = idx[j],
        key = x + "," + y;
      edges.set(key, (edges.get(key) || 0) + 1);
    }
  }
  for (const [key, count] of edges)
    if (count >= opts.support) {
      const [a, b] = key.split(",").map(Number);
      if (!variants[a].overlap && !variants[b].overlap)
        parent[root(b)] = root(a);
    }
  const components = new Map();
  variants.forEach((v, i) => {
    const k = root(i);
    if (!components.has(k)) components.set(k, []);
    components.get(k).push(i);
  });
  const blocks = [];
  const assignments = new Map();
  for (const indices of components.values()) {
    const patterns = new Map();
    for (let k = 0; k < fs.length; k++) {
      const pattern = indices.map((i) => vectors[k][i]);
      if (!pattern.some((a) => a !== "?")) continue;
      const key = JSON.stringify(pattern);
      if (!patterns.has(key)) patterns.set(key, { pattern, members: [] });
      patterns.get(key).members.push(k);
    }
    // Exact observed patterns remain separate: missing values never invent a full haplotype.
    const seeds = [...patterns.values()]
      .filter((p) => p.members.length >= opts.support)
      .sort(
        (a, b) =>
          b.pattern.filter((x) => x !== "?").length -
            a.pattern.filter((x) => x !== "?").length ||
          b.members.length - a.members.length,
      );
    const groups = [];
    for (const seed of seeds) {
      if (
        !groups.some((g) =>
          seed.pattern.every((a, i) => a === "?" || a === g.pattern[i]),
        )
      )
        groups.push({ ...seed, members: [], compatible: 0 });
    }
    if (groups.length > 32)
      throw Error(
        "Mixture exceeds 32 supported patterns; increase support or narrow the interval",
      );
    let ambiguous = 0,
      conflicting = 0;
    for (let k = 0; k < fs.length; k++) {
      const a = indices.map((i) => vectors[k][i]);
      if (!a.some((x) => x !== "?")) continue;
      const match = groups.filter((g) =>
        a.every(
          (v, i) => v === "?" || g.pattern[i] === "?" || v === g.pattern[i],
        ),
      );
      if (
        match.length === 1 &&
        a.some((v, i) => v !== "?" && match[0].pattern[i] !== "?")
      ) {
        match[0].members.push(k);
        const list = assignments.get(k) || [];
        list.push(`B${blocks.length + 1} H${groups.indexOf(match[0]) + 1}`);
        assignments.set(k, list);
      } else if (match.length > 1) {
        ambiguous++;
        match.forEach((g) => g.compatible++);
      } else conflicting++;
    }
    blocks.push({
      id: blocks.length + 1,
      indices,
      start: variants[indices[0]].pos,
      end:
        variants[indices.at(-1)].pos +
        Math.max(1, variants[indices.at(-1)].ref.length),
      groups: groups.map((g, i) => ({
        ...g,
        id: i + 1,
        count: g.members.length,
      })),
      ambiguous,
      conflicting,
    });
  }
  const coverage = new Uint32Array(reference.length);
  for (const f of fs)
    for (const p of f.bases.keys())
      if (p >= origin && p < origin + reference.length) coverage[p - origin]++;
  const events = [...counts.values()]
    .filter((v) => v.ref.length !== v.alt.length && v.support >= opts.support)
    .map((v) => {
      let type = v.ref.length ? "Deletion" : "Insertion";
      let source = null;
      if (!v.ref.length) {
        const local = reference.slice(
          Math.max(0, v.pos - origin - v.alt.length),
          v.pos - origin,
        );
        if (local === v.alt) {
          type = "Candidate tandem duplication";
          source = [v.pos - v.alt.length, v.pos];
        } else if (
          reference.slice(v.pos - origin, v.pos - origin + v.alt.length) ===
          v.alt
        ) {
          type = "Candidate tandem duplication";
          source = [v.pos, v.pos + v.alt.length];
        }
      }
      return { ...v, type, source };
    });
  const clipMap = new Map();
  for (const f of fs) {
    const seen = new Set();
    for (const r of f.reads)
      for (const c of r.projection.clips) {
        if (c.seq.length < 8 || c.seq.length > 1000) continue;
        const key = JSON.stringify([c.pos, c.side, c.seq]);
        if (seen.has(key)) continue;
        seen.add(key);
        const v = clipMap.get(key) || { ...c, support: 0 };
        v.support++;
        clipMap.set(key, v);
      }
  }
  return {
    variants,
    blocks,
    coverage: [...coverage],
    excluded,
    fragmentCount: fs.length,
    events,
    clips: [...clipMap.values()]
      .filter((c) => c.support >= opts.support)
      .sort((a, b) => b.support - a.support)
      .slice(0, 32),
    vectors,
    fragments: fs.map((f, k) => ({
      id: k,
      reads: f.reads.map((r) => r.id),
      groups: assignments.get(k) || [],
    })),
    reads: reads.map((r) => {
      const { projection, ...rest } = r;
      return rest;
    }),
    edges: [...edges],
  };
}
// Bounded affine-gap global alignment. Used on anchored local windows, never a whole gene.
export function align(a, b, limit = 2e6) {
  const n = a.length,
    m = b.length;
  if ((n + 1) * (m + 1) > limit)
    throw Error(
      "Alignment cell budget exceeded; select a smaller event window",
    );
  const size = (n + 1) * (m + 1),
    M = new Float32Array(size),
    I = new Float32Array(size),
    D = new Float32Array(size),
    tm = new Uint8Array(size),
    ti = new Uint8Array(size),
    td = new Uint8Array(size);
  M.fill(-1e9);
  I.fill(-1e9);
  D.fill(-1e9);
  M[0] = 0;
  for (let i = 1; i <= n; i++) D[i * (m + 1)] = -5 - i;
  for (let j = 1; j <= m; j++) I[j] = -5 - j;
  for (let i = 1; i <= n; i++)
    for (let j = 1; j <= m; j++) {
      const k = i * (m + 1) + j,
        l = k - 1,
        u = k - m - 1,
        z = u - 1,
        prev = [M[z], I[z], D[z]];
      tm[k] = prev.indexOf(Math.max(...prev));
      M[k] = Math.max(...prev) + (a[i - 1] === b[j - 1] ? 2 : -3);
      const ip = [M[l] - 6, I[l] - 1, D[l] - 6];
      ti[k] = ip.indexOf(Math.max(...ip));
      I[k] = Math.max(...ip);
      const dp = [M[u] - 6, I[u] - 6, D[u] - 1];
      td[k] = dp.indexOf(Math.max(...dp));
      D[k] = Math.max(...dp);
    }
  let i = n,
    j = m,
    k = n * (m + 1) + m,
    state = [M[k], I[k], D[k]].indexOf(Math.max(M[k], I[k], D[k])),
    x = "",
    y = "";
  const score = Math.max(M[k], I[k], D[k]);
  while (i || j) {
    k = i * (m + 1) + j;
    if (!i) state = 1;
    if (!j) state = 2;
    if (state === 0) {
      x = a[--i] + x;
      y = b[--j] + y;
      state = tm[k];
    } else if (state === 1) {
      x = "-" + x;
      y = b[--j] + y;
      state = ti[k];
    } else {
      x = a[--i] + x;
      y = "-" + y;
      state = td[k];
    }
  }
  return { a: x, b: y, score };
}
export function eventPath(event, reference, origin, flank = 60) {
  const start = Math.max(origin, event.pos - flank),
    end = Math.min(
      origin + reference.length,
      event.pos + event.ref.length + flank,
    );
  const ref = reference.slice(start - origin, end - origin),
    offset = event.pos - start,
    alt =
      ref.slice(0, offset) + event.alt + ref.slice(offset + event.ref.length);
  const map = [];
  for (let i = start; i < event.pos; i++) map.push({ position: i });
  for (let i = 0; i < event.alt.length; i++)
    map.push(
      event.ref.length
        ? { position: event.pos + i }
        : { anchor: event.pos, offset: i },
    );
  for (let i = event.pos + event.ref.length; i < end; i++)
    map.push({ position: i });
  return { start, end, ref, alt, map, alignment: align(ref, alt) };
}

export function cropRead(read, start, end) {
  let p = read.start,
    q = 0,
    out = "";
  for (const [n, op] of cigarOps(read.cigar)) {
    if ("M=X".includes(op)) {
      const a = Math.max(0, start - p),
        b = Math.min(n, end - p);
      if (a < b) out += read.seq.slice(q + a, q + b);
      p += n;
      q += n;
    } else if (op === "I") {
      if (p >= start && p < end) out += read.seq.slice(q, q + n);
      q += n;
    } else if (op === "S") q += n;
    else if (op === "D" || op === "N") p += n;
  }
  return out;
}
