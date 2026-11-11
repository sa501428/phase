import { recordEvidence } from "./model.js";
import { normalize } from "../src/analysis.js";

const popcount = (x) => {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  return (((x + (x >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
};
const intersects = (a, b) => {
  let n = 0;
  for (let i = 0; i < a.length; i++) n += popcount(a[i] & b[i]);
  return n;
};
const put = (bits, i) => {
  bits[i >>> 5] |= 1 << (i & 31);
};
const short = (s) => (s.length > 10 ? `${s.slice(0, 7)}…(${s.length} bp)` : s);
export const variantLabel = (v) =>
  v.kind === "insertion"
    ? `${v.pos + 1} ins ${short(v.alt)}`
    : v.kind === "deletion"
      ? `${v.pos + 1} del ${v.ref.length} bp`
      : `${v.pos + 1} ${v.ref}>${v.alt}`;

// A fragment is the unit of observation. Both mates may bridge different sites,
// but conflicting mates and absent coverage never become reference evidence.
export function buildVariantGraph(model, { maxVariants = 180 } = {}) {
  const ids = model.activeIDs(model.fragments.map((f) => f.id)),
    fragments = [],
    candidates = new Map();
  const candidate = (v) => {
    if (!candidates.has(v.id))
      candidates.set(v.id, { ...v, anchors: new Set() });
    return candidates.get(v.id);
  };
  const hasIndels = model.reads.some((r) =>
    (r.runs || []).some(
      (run) =>
        ["I", "D"].includes(run.op) &&
        run.refStart >= model.region.start &&
        run.refStart < model.region.end,
    ),
  );
  for (const id of ids) {
    const evidence = model.evidence(id),
      records = [];
    for (const [p, values] of evidence.bases) {
      const ref = model.reference[p - model.region.start];
      if (!/^[ACGT]$/.test(ref) || values.size !== 1) continue;
      const alt = [...values][0];
      if (/^[ACGT]$/.test(alt) && alt !== ref)
        candidate({ id: `${p}:${ref}>${alt}`, pos: p, ref, alt, kind: "SNV" });
    }
    if (!hasIndels) {
      fragments.push({ id, evidence, records });
      continue;
    }
    for (const rid of model.fragments[id].readIds) {
      const r = model.readById.get(rid),
        e = recordEvidence(r, model.region, model.settings),
        events = [];
      for (const run of r.runs || []) {
        const { op, refStart: p, queryStart: q, length: n } = run;
        if (
          !["I", "D"].includes(op) ||
          p < model.region.start ||
          p >= model.region.end ||
          (op === "D" && p + n > model.region.end)
        )
          continue;
        const ref =
            op === "D"
              ? model.reference.slice(
                  p - model.region.start,
                  p - model.region.start + n,
                )
              : "",
          alt = op === "I" ? r.seq?.slice(q, q + n) : "";
        if (!/^[ACGT]*$/.test(ref + (alt || "")) || (op === "I" && !alt))
          continue;
        const v = {
          ...normalize(p, ref, alt, model.reference, model.region.start),
          kind: op === "I" ? "insertion" : "deletion",
        };
        const valid =
          op === "I"
            ? e.insertions.get(p) === alt
            : Array.from(
                { length: n },
                (_, j) => e.bases.get(p + j) === "-",
              ).every(Boolean);
        events.push({ ...v, valid, rawPos: p });
        if (valid) candidate(v).anchors.add(p);
      }
      records.push({ e, events, runs: r.runs || [] });
    }
    fragments.push({ id, evidence, records });
  }
  if (candidates.size > maxVariants)
    return {
      nodes: [],
      edges: [],
      fragments: ids.length,
      error: `${candidates.size} visible alleles exceed the graph limit of ${maxVariants}. Narrow the loaded region or raise the VAF cutoff. All reads remain available.`,
    };
  const raw = [...candidates.values()].sort(
      (a, b) => a.pos - b.pos || a.id.localeCompare(b.id),
    ),
    nodes = [],
    words = Math.ceil(ids.length / 32);
  for (const v of raw) {
    const altBits = new Uint32Array(words),
      refBits = new Uint32Array(words),
      carriers = [];
    let altCount = 0,
      total = 0,
      other = 0,
      conflicting = 0;
    for (let i = 0; i < fragments.length; i++) {
      const f = fragments[i];
      let state = "?";
      if (v.kind === "SNV") {
        const calls = f.evidence.bases.get(v.pos);
        if (calls?.size > 1) state = "X";
        else if (calls?.size === 1) {
          const b = [...calls][0];
          state = b === v.alt ? "A" : b === v.ref ? "R" : "O";
        }
      } else {
        const votes = new Set();
        for (const r of f.records) {
          const exact = r.events.filter((e) => e.id === v.id);
          if (exact.length) {
            if (exact.some((e) => e.valid)) votes.add("A");
            continue;
          }
          // An alternative or filtered indel at this boundary is not REF.
          const overlap = r.events.some(
            (e) =>
              e.pos === v.pos ||
              (e.kind === "deletion" &&
                e.pos < v.pos + Math.max(1, v.ref.length) &&
                e.pos + e.ref.length > v.pos),
          );
          if (overlap) continue;
          const spanned = [v.pos, ...v.anchors].some((p) => {
            const end = p + v.ref.length;
            return (
              r.runs.some(
                (run) =>
                  "M=X".includes(run.op) &&
                  run.refStart < p &&
                  run.refStart + run.length > end,
              ) &&
              r.e.bases.has(p - 1) &&
              r.e.bases.has(end)
            );
          });
          if (spanned) votes.add("R");
        }
        state = votes.size > 1 ? "X" : votes.size === 1 ? [...votes][0] : "?";
      }
      if (state === "X") conflicting++;
      else if (state === "O") {
        other++;
        total++;
      } else if (state === "A" || state === "R") {
        total++;
        put(state === "A" ? altBits : refBits, i);
        if (state === "A") {
          altCount++;
          carriers.push(f.id);
        }
      }
    }
    const vaf = total ? altCount / total : 0;
    if (!altCount || vaf < (model.settings.minVaf ?? 0.02)) continue;
    nodes.push({
      id: v.id,
      index: nodes.length,
      pos: v.pos,
      ref: v.ref,
      alt: v.alt,
      kind: v.kind,
      label: variantLabel(v),
      altCount,
      total,
      vaf,
      other,
      conflicting,
      unknown: fragments.length - total - conflicting,
      carriers,
      altBits,
      refBits,
    });
  }
  const edges = [];
  for (let a = 0; a < nodes.length; a++)
    for (let b = a + 1; b < nodes.length; b++) {
      const x = nodes[a],
        y = nodes[b],
        both = intersects(x.altBits, y.altBits),
        leftOnly = intersects(x.altBits, y.refBits),
        rightOnly = intersects(x.refBits, y.altBits),
        neither = intersects(x.refBits, y.refBits),
        callable = both + leftOnly + rightOnly + neither;
      if (!both) continue;
      edges.push({
        id: `${a}:${b}`,
        a,
        b,
        both,
        leftOnly,
        rightOnly,
        neither,
        callable,
        unknown: fragments.length - callable,
        frequency: both / callable,
        linkage: both / (both + leftOnly + rightOnly),
      });
    }
  for (const n of nodes) {
    delete n.altBits;
    delete n.refBits;
  }
  return { nodes, edges, fragments: ids.length };
}

// Contrast is relative to the strongest displayed connection, with an explicit
// legend. The same metric controls both thickness and opacity.
export function edgeAppearance(value, maximum) {
  const t = maximum ? Math.max(0, Math.min(1, value / maximum)) : 0;
  return {
    width: 0.65 + 16 * Math.pow(t, 1.55),
    opacity: 0.13 + 0.82 * Math.sqrt(t),
  };
}
export function layoutVariantGraph(nodes, edges) {
  const points = nodes.map((n, i) => ({
    x: 500 + 220 * Math.cos(i * 2.399963),
    y: 210 + 150 * Math.sin(i * 2.399963),
    vx: 0,
    vy: 0,
  }));
  const adjacency = nodes.map(() => []);
  for (const e of edges) {
    adjacency[e.a].push(e.b);
    adjacency[e.b].push(e.a);
  }
  const seen = new Set(),
    components = [];
  for (let i = 0; i < nodes.length; i++)
    if (!seen.has(i)) {
      const c = [i];
      seen.add(i);
      for (let k = 0; k < c.length; k++)
        for (const j of adjacency[c[k]])
          if (!seen.has(j)) {
            seen.add(j);
            c.push(j);
          }
      components.push(c);
    }
  // Separate disconnected coverage blocks, including isolated variants.
  components.forEach((c, k) => {
    const center = {
      x:
        500 +
        (components.length > 1
          ? 300 * Math.cos((k * 2 * Math.PI) / components.length)
          : 0),
      y:
        210 +
        (components.length > 1
          ? 110 * Math.sin((k * 2 * Math.PI) / components.length)
          : 0),
    };
    c.forEach((i) => {
      points[i].cx = center.x;
      points[i].cy = center.y;
    });
  });
  for (let step = 0; step < 180; step++) {
    const forces = points.map((p) => ({
      x: (p.cx - p.x) * 0.012,
      y: (p.cy - p.y) * 0.018,
    }));
    for (let a = 0; a < points.length; a++)
      for (let b = a + 1; b < points.length; b++) {
        const dx = points[b].x - points[a].x,
          dy = points[b].y - points[a].y,
          d2 = Math.max(160, dx * dx + dy * dy),
          f = 2600 / d2;
        forces[a].x -= dx * f;
        forces[a].y -= dy * f;
        forces[b].x += dx * f;
        forces[b].y += dy * f;
      }
    for (const e of edges) {
      const p = points[e.a],
        q = points[e.b],
        dx = q.x - p.x,
        dy = q.y - p.y,
        dist = Math.hypot(dx, dy) || 1;
      const f =
        ((dist - (65 + 105 * (1 - e.linkage))) * (0.008 + 0.025 * e.linkage)) /
        dist;
      forces[e.a].x += dx * f;
      forces[e.a].y += dy * f;
      forces[e.b].x -= dx * f;
      forces[e.b].y -= dy * f;
    }
    points.forEach((p, i) => {
      p.vx = (p.vx + forces[i].x) * 0.66;
      p.vy = (p.vy + forces[i].y) * 0.66;
      p.x += p.vx;
      p.y += p.vy;
    });
  }
  return points;
}
