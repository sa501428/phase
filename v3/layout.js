// IGV-style interval partitioning: a pair occupies its complete template span,
// including the connector and visible soft clips. Nonoverlapping templates can
// reuse a lane; overlapping templates always remain separate selectable items.
class MinHeap {
  constructor() {
    this.values = [];
  }
  push(value) {
    const a = this.values;
    a.push(value);
    let i = a.length - 1;
    while (i) {
      const p = (i - 1) >> 1;
      if (
        a[p].end < value.end ||
        (a[p].end === value.end && a[p].lane <= value.lane)
      )
        break;
      a[i] = a[p];
      i = p;
    }
    a[i] = value;
  }
  pop() {
    const a = this.values,
      first = a[0],
      last = a.pop();
    if (a.length) {
      let i = 0;
      while (i * 2 + 1 < a.length) {
        let j = i * 2 + 1;
        if (
          j + 1 < a.length &&
          (a[j + 1].end < a[j].end ||
            (a[j + 1].end === a[j].end && a[j + 1].lane < a[j].lane))
        )
          j++;
        if (
          last.end < a[j].end ||
          (last.end === a[j].end && last.lane <= a[j].lane)
        )
          break;
        a[i] = a[j];
        i = j;
      }
      a[i] = last;
    }
    return first;
  }
  peek() {
    return this.values[0];
  }
}
export function fragmentSpan(model, f) {
  let start = Infinity,
    end = -Infinity;
  for (const id of f.readIds) {
    const r = model.readById.get(id);
    if (r.chr && r.chr !== model.region.chr) continue;
    if (r.alignmentUsable) {
      start = Math.min(start, r.start);
      end = Math.max(end, r.end);
      for (const run of r.runs)
        if (run.op === "S") {
          const s =
            run.queryStart === 0 ? run.refStart - run.length : run.refStart;
          start = Math.min(start, s);
          end = Math.max(end, s + run.length);
        }
    } else if (Number.isFinite(r.start) && r.start >= 0) {
      start = Math.min(start, r.start);
      end = Math.max(end, r.start + 1);
    }
  }
  return { start, end };
}
export function packPileup(
  model,
  items,
  {
    view,
    width,
    left = 112,
    readHeight = 18,
    sort = "position",
    site = null,
    minHeight = 570,
  },
) {
  const stride = readHeight + 9,
    placements = [],
    byItem = new Map(),
    sections = [];
  let top = 8,
    laneCount = 0;
  const scale = (width - left - 12) / (view.end - view.start),
    x = (p) => left + (p - view.start) * scale;
  const place = (index, y, height, span, lane, header = false) => {
    const p = {
      index,
      key: items[index].key,
      y,
      height,
      start: span.start,
      end: span.end,
      x1: Math.max(left, x(span.start)),
      x2: Math.min(width - 12, Math.max(x(span.end), x(span.start) + 8)),
      lane,
      header,
    };
    placements.push(p);
    byItem.set(index, p);
  };
  function pack(indices, label) {
    const visible = indices
      .map((index) => ({
        index,
        span: fragmentSpan(model, items[index].fragment),
      }))
      .filter((p) => p.span.end > view.start && p.span.start < view.end);
    if (!visible.length) return;
    const sectionTop = top;
    const bands = new Map();
    for (const p of visible) {
      let band = "";
      if (sort === "allele" && site !== null) {
        const values = model
          .evidence(items[p.index].fragment.id)
          .bases.get(site);
        band =
          values?.size === 1 ? [...values][0] : values?.size > 1 ? "!" : "?";
      }
      if (!bands.has(band)) bands.set(band, []);
      bands.get(band).push(p);
    }
    for (const band of [...bands.keys()].sort()) {
      const entries = bands
        .get(band)
        .sort(
          (a, b) =>
            a.span.start - b.span.start ||
            (sort === "name"
              ? items[a.index].fragment.name.localeCompare(
                  items[b.index].fragment.name,
                )
              : a.index - b.index),
        );
      const heap = new MinHeap();
      let lanes = 0;
      for (const entry of entries) {
        const start = Math.max(left, x(entry.span.start)) - 5,
          end =
            Math.min(
              width - 12,
              Math.max(x(entry.span.end), x(entry.span.start) + 8),
            ) + 5;
        let lane;
        if (heap.peek()?.end <= start) lane = heap.pop().lane;
        else lane = lanes++;
        heap.push({ end, lane });
        place(entry.index, top + lane * stride, readHeight, entry.span, lane);
      }
      top += lanes * stride;
      laneCount += lanes;
    }
    sections.push({ label, y: sectionTop, end: top });
    top += 9;
  }
  let pending = [],
    pendingLabel = "Read pairs";
  const flush = () => {
    pack(pending, pendingLabel);
    pending = [];
  };
  for (let i = 0; i < items.length; i++) {
    const item = items[i];
    if (item.kind === "group") {
      flush();
      place(i, top, Math.max(30, readHeight + 12), view, -1, true);
      top += Math.max(30, readHeight + 12) + 5;
      pendingLabel = "↳ Members";
    } else {
      if (item.kind === "fragment" && pendingLabel === "↳ Members") {
        flush();
        pendingLabel = "Read pairs";
      }
      pending.push(i);
    }
  }
  flush();
  placements.sort((a, b) => a.y - b.y || a.x1 - b.x1 || a.index - b.index);
  return {
    placements,
    byItem,
    sections,
    order: placements.map((p) => p.index),
    height: Math.max(minHeight, top),
    laneCount,
    readHeight,
    stride,
  };
}
export function pileupHit(layout, x, y) {
  return (
    layout.placements.find(
      (p) =>
        y >= p.y - 2 &&
        y <= p.y + p.height + 2 &&
        x >= p.x1 - 4 &&
        x <= p.x2 + 4,
    ) || null
  );
}
