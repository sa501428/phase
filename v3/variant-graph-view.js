import {
  buildVariantGraph,
  layoutVariantGraph,
  edgeAppearance,
} from "./variant-graph.js";
const ns = "http://www.w3.org/2000/svg";
const svg = (name, attrs = {}) => {
  const e = document.createElementNS(ns, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
};
const text = (tag, value, cls) => {
  const e = document.createElement(tag);
  e.textContent = value;
  if (cls) e.className = cls;
  return e;
};
const pct = (n) => (n * 100).toFixed(1) + "%";
const fullLabel = (n, chr) =>
  n.kind === "insertion"
    ? `${chr}:${n.pos + 1} ins ${n.alt}`
    : n.kind === "deletion"
      ? `${chr}:${n.pos + 1} del ${n.ref}`
      : `${chr}:${n.pos + 1}${n.ref}>${n.alt}`;
export class VariantGraphView {
  constructor(panel, callbacks) {
    this.panel = panel;
    this.callbacks = callbacks;
    this.model = null;
    this.data = null;
    this.dirty = true;
    this.selection = null;
    this.$ = (id) => panel.querySelector("#" + id);
    this.plot = this.$("variantNetwork");
    panel.addEventListener("toggle", () => {
      if (panel.open) this.ensure();
    });
    for (const id of ["graphWeight", "graphSupport"])
      this.$(id).addEventListener("change", () => this.draw());
    this.$("graphFit").onclick = () => this.fit();
    this.$("graphIn").onclick = () => this.zoom(0.75);
    this.$("graphOut").onclick = () => this.zoom(1 / 0.75);
    this.plot.addEventListener(
      "wheel",
      (e) => {
        if (!this.data?.nodes.length) return;
        e.preventDefault();
        this.zoom(e.deltaY > 0 ? 1.15 : 1 / 1.15);
      },
      { passive: false },
    );
    this.plot.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 || !this.data?.nodes.length) return;
      const node = e.target.closest("[data-node]");
      if (e.target.closest("[data-edge]")) return;
      this.drag = {
        x: e.clientX,
        y: e.clientY,
        box: [...this.box],
        node: node ? +node.dataset.node : null,
        moved: false,
      };
      (node || this.plot).setPointerCapture(e.pointerId);
    });
    this.plot.addEventListener("pointermove", (e) => {
      if (!this.drag) return;
      const b = this.plot.getBoundingClientRect(),
        scale = Math.min(b.width / this.box[2], b.height / this.box[3]),
        dx = (e.clientX - this.drag.x) / scale,
        dy = (e.clientY - this.drag.y) / scale;
      if (Math.hypot(e.clientX - this.drag.x, e.clientY - this.drag.y) > 3)
        this.drag.moved = true;
      if (this.drag.node !== null) {
        const p = this.points[this.drag.node];
        if (!this.drag.point) this.drag.point = { x: p.x, y: p.y };
        p.x = this.drag.point.x + dx;
        p.y = this.drag.point.y + dy;
        this.position();
      } else {
        this.box = [
          this.drag.box[0] - dx,
          this.drag.box[1] - dy,
          this.drag.box[2],
          this.drag.box[3],
        ];
        this.viewBox();
      }
    });
    this.plot.addEventListener("pointerup", () => {
      this.suppressClick = this.drag?.moved;
      this.drag = null;
    });
    this.plot.addEventListener("pointercancel", () => {
      this.drag = null;
    });
    this.plot.addEventListener("click", (e) => {
      if (this.suppressClick) {
        this.suppressClick = false;
        return;
      }
      const n = e.target.closest("[data-node]"),
        edge = e.target.closest("[data-edge]");
      if (n) this.choose({ kind: "node", index: +n.dataset.node });
      else if (edge) this.choose({ kind: "edge", index: +edge.dataset.edge });
      else this.choose(null);
    });
    this.plot.addEventListener("keydown", (e) => {
      const n = e.target.closest("[data-node]"),
        edge = e.target.closest("[data-edge]");
      if ((e.key === "Enter" || e.key === " ") && (n || edge)) {
        e.preventDefault();
        this.choose(
          n
            ? { kind: "node", index: +n.dataset.node }
            : { kind: "edge", index: +edge.dataset.edge },
        );
      } else if (e.key === "Escape") this.choose(null);
    });
    new ResizeObserver(() => {
      if (this.panel.open && this.data) this.viewBox();
    }).observe(this.plot);
    this.update(null);
  }
  update(model) {
    this.drag = null;
    this.suppressClick = false;
    this.model = model;
    this.dirty = true;
    this.selection = null;
    if (!model || this.panel.open) this.ensure();
    else
      this.$("graphSummary").textContent = model
        ? "Expand to view variant co-occurrence"
        : "Load a region to view variant connections";
  }
  ensure() {
    if (!this.dirty) return;
    this.dirty = false;
    this.data = this.model
      ? buildVariantGraph(this.model)
      : { nodes: [], edges: [], fragments: 0 };
    this.points = layoutVariantGraph(this.data.nodes, this.data.edges);
    this.draw();
    this.fit();
  }
  draw() {
    if (!this.data) return;
    const { nodes, edges } = this.data,
      mode = this.$("graphWeight").value,
      min = Math.max(1, Math.floor(+this.$("graphSupport").value) || 1);
    this.metric = mode === "linkage" ? "linkage" : "frequency";
    this.visible = edges
      .map((e, index) => ({ ...e, index }))
      .filter((e) => e.both >= min);
    this.maximum = Math.max(0, ...this.visible.map((e) => e[this.metric]));
    const minimum = this.visible.length
      ? Math.min(...this.visible.map((e) => e[this.metric]))
      : 0;
    this.plot.replaceChildren();
    this.nodeElements = [];
    this.edgeElements = [];
    this.plot.dataset.nodes = nodes.length;
    this.plot.dataset.edges = this.visible.length;
    this.plot.dataset.minVaf = this.model?.settings.minVaf ?? "";
    this.$("graphSummary").textContent = this.data.error
      ? "Graph limit reached"
      : `${nodes.length} variants · ${this.visible.length} connections · ${this.data.fragments} included pairs`;
    this.$("graphMessage").textContent =
      this.data.error ||
      (!this.model
        ? "Load a region to view variant connections."
        : !nodes.length
          ? "No callable SNPs or indels meet the current VAF cutoff."
          : "");
    this.$("graphMessage").hidden = !this.$("graphMessage").textContent;
    this.$("graphScale").textContent = this.visible.length
      ? `Thin ${pct(minimum)} → thick ${pct(this.maximum)} (relative to this graph)`
      : "No measured connections";
    this.$("graphDefinition").textContent =
      this.metric === "frequency"
        ? "Line weight = both ALT ÷ pairs callable as REF/ALT at both sites."
        : "Line weight = both ALT ÷ (both ALT + left-only ALT + right-only ALT).";
    for (const e of this.visible) {
      const style = edgeAppearance(e[this.metric], this.maximum),
        g = svg("g", {
          "data-edge": e.index,
          tabindex: "0",
          role: "button",
          "aria-label": `Connection V${e.a + 1} to V${e.b + 1}: ${e.both}/${e.callable} both ALT`,
        }),
        title = svg("title");
      title.textContent = `V${e.a + 1} ↔ V${e.b + 1}: both ALT ${e.both}/${e.callable} (${pct(e.frequency)}); ALT linkage ${pct(e.linkage)}`;
      const hit = svg("line", {
          stroke: "transparent",
          "vector-effect": "non-scaling-stroke",
          "stroke-width": Math.max(10, style.width),
          "pointer-events": "stroke",
        }),
        line = svg("line", {
          stroke: "#207d80",
          "vector-effect": "non-scaling-stroke",
          "stroke-width": style.width,
          "stroke-opacity": style.opacity,
          "stroke-linecap": "round",
          "pointer-events": "none",
          "data-weight": e[this.metric],
          "data-shared": e.both,
        });
      g.append(title, hit, line);
      this.plot.append(g);
      this.edgeElements.push({ e, g, line, hit });
    }
    for (const n of nodes) {
      const g = svg("g", {
          "data-node": n.index,
          tabindex: "0",
          role: "button",
          "aria-label": `V${n.index + 1}: ${fullLabel(n, this.model.region.chr)}, VAF ${pct(n.vaf)}`,
        }),
        title = svg("title");
      title.textContent = `V${n.index + 1} · ${fullLabel(n, this.model.region.chr)}\nALT ${n.altCount}/${n.total} callable pairs · VAF ${pct(n.vaf)}`;
      const circle = svg("circle", {
          r: 10 + 16 * Math.sqrt(n.vaf),
          fill: { SNV: "#d8edf6", insertion: "#e9ddf7", deletion: "#fae4cd" }[
            n.kind
          ],
          stroke: { SNV: "#37758c", insertion: "#8b5daa", deletion: "#b47b39" }[
            n.kind
          ],
          "stroke-width": 2,
        }),
        label = svg("text", {
          "text-anchor": "middle",
          dy: 4,
          class: "graph-node-number",
        });
      label.textContent = "V" + (n.index + 1);
      const caption = svg("text", {
        "text-anchor": "middle",
        y: 34,
        class: "graph-node-label",
      });
      caption.textContent = n.label;
      g.append(title, circle, label, caption);
      this.plot.append(g);
      this.nodeElements.push(g);
    }
    this.position();
    this.choose(null);
  }
  position() {
    for (const { e, line, hit } of this.edgeElements)
      for (const el of [line, hit])
        for (const [attr, value] of Object.entries({
          x1: this.points[e.a].x,
          y1: this.points[e.a].y,
          x2: this.points[e.b].x,
          y2: this.points[e.b].y,
        }))
          el.setAttribute(attr, value);
    this.positionNodes();
  }
  positionNodes() {
    if (!this.box) return;
    const bounds = this.plot.getBoundingClientRect();
    const scale =
      Math.min(bounds.width / this.box[2], bounds.height / this.box[3]) || 1;
    this.nodeElements.forEach((g, i) =>
      g.setAttribute(
        "transform",
        `translate(${this.points[i].x} ${this.points[i].y}) scale(${1 / scale})`,
      ),
    );
  }
  fit() {
    if (!this.points?.length) {
      this.box = [0, 0, 1000, 350];
      this.viewBox();
      return;
    }
    const xs = this.points.map((p) => p.x),
      ys = this.points.map((p) => p.y),
      bounds = this.plot.getBoundingClientRect();
    const scale = Math.min(
      Math.max(60, bounds.width - 180) /
        Math.max(160, Math.max(...xs) - Math.min(...xs)),
      Math.max(50, bounds.height - 110) /
        Math.max(100, Math.max(...ys) - Math.min(...ys)),
      1.8,
    );
    const w = Math.max(100, bounds.width) / scale,
      h = Math.max(100, bounds.height) / scale;
    this.box = [
      (Math.min(...xs) + Math.max(...xs) - w) / 2,
      (Math.min(...ys) + Math.max(...ys) - h) / 2 + 8 / scale,
      w,
      h,
    ];
    this.viewBox();
  }
  viewBox() {
    if (this.box) {
      this.plot.setAttribute("viewBox", this.box.join(" "));
      this.positionNodes();
    }
  }
  zoom(factor) {
    if (!this.box) return;
    const [x, y, w, h] = this.box,
      nw = Math.max(100, Math.min(6000, w * factor)),
      nh = (h * nw) / w;
    this.box = [x + (w - nw) / 2, y + (h - nh) / 2, nw, nh];
    this.viewBox();
  }
  choose(selection) {
    this.selection = selection;
    const detail = this.$("graphDetail");
    detail.replaceChildren();
    let linked = new Set(),
      carriers = [],
      focus = null;
    if (selection?.kind === "node") {
      const n = this.data.nodes[selection.index];
      linked.add(n.index);
      carriers = n.carriers;
      focus = n;
      for (const e of this.visible)
        if (e.a === n.index || e.b === n.index) {
          linked.add(e.a);
          linked.add(e.b);
        }
      detail.append(
        text("strong", `V${n.index + 1} · ${n.label}`),
        text(
          "span",
          `ALT ${n.altCount}/${n.total} callable pairs (${pct(n.vaf)} VAF). ${n.conflicting} conflicting pairs; ${n.other} other alleles; ${n.unknown} unknown.`,
        ),
      );
      const carrierSet = new Set(carriers);
      const groups = this.model.groups
        .map((g) => ({
          name: g.name,
          n: g.members.filter((id) => carrierSet.has(id)).length,
        }))
        .filter((g) => g.n);
      if (groups.length)
        detail.append(
          text(
            "span",
            "Constructed groups: " +
              groups.map((g) => `${g.name} (${g.n})`).join(" · "),
          ),
        );
    } else if (selection?.kind === "edge") {
      const e = this.data.edges[selection.index],
        a = this.data.nodes[e.a],
        b = this.data.nodes[e.b],
        bs = new Set(b.carriers);
      linked = new Set([e.a, e.b]);
      carriers = a.carriers.filter((id) => bs.has(id));
      focus = a;
      detail.append(
        text("strong", `V${e.a + 1} ↔ V${e.b + 1}`),
        text(
          "span",
          `Both ALT ${e.both}/${e.callable} (${pct(e.frequency)} joint frequency; ${pct(e.linkage)} ALT linkage). A-only ${e.leftOnly}; B-only ${e.rightOnly}; neither ${e.neither}; unknown/other/conflicting ${e.unknown}.`,
        ),
      );
    } else
      detail.append(
        text(
          "span",
          (this.data.nodes.length && !this.visible.length
            ? "No shared ALT pairs meet the connection setting. "
            : "") +
            "Click a node or line for counts. Drag nodes to arrange; drag the background to pan; wheel to zoom.",
          "muted",
        ),
      );
    if (selection) {
      const select = text(
        "button",
        `Select ${carriers.length} supporting pairs`,
      );
      select.type = "button";
      select.onclick = () => this.callbacks.selectPairs(carriers);
      const locate = text("button", "Show site in pileup");
      locate.type = "button";
      locate.onclick = () => this.callbacks.focusVariant(focus);
      detail.append(select, locate);
    }
    for (const [i, g] of this.nodeElements.entries()) {
      g.classList.toggle("graph-muted", !!selection && !linked.has(i));
      g.classList.toggle(
        "graph-chosen",
        selection?.kind === "node" && selection.index === i,
      );
    }
    for (const { e, g } of this.edgeElements) {
      const active =
        !selection ||
        (selection.kind === "node"
          ? e.a === selection.index || e.b === selection.index
          : e.index === selection.index);
      g.classList.toggle("graph-muted", !active);
      g.classList.toggle("graph-chosen", !!selection && active);
    }
  }
}
