import { watchVaf } from "../src/live-vaf.js";
import { hiddenBase, hiddenEvent } from "../src/vaf-filter.js";
import { ungzip } from "pako";
import hg19 from "../references/hg19.js";
import hg38 from "../references/hg38.js";
import { alignmentRuns, hasSequence } from "../src/alignment-geometry.js";
import { reverseComplement } from "./assembly.js";

const $ = (id) => document.getElementById(id),
  make = (tag, text, cls) => {
    const e = document.createElement(tag);
    if (text !== undefined) e.textContent = text;
    if (cls) e.className = cls;
    return e;
  };
let worker,
  workerURL,
  timer,
  opened = false,
  busy = false,
  result = null,
  contigs = [],
  catalog = null,
  catalogBuild,
  genes = [],
  selected = null,
  selection = null,
  query = null;
let positions = [],
  camera = { x: 0, y: 0, scale: 1 },
  highlight = new Set(),
  queried = new Set(),
  readFragments = [],
  readView = null,
  readsOpen = false;
const graph = $("graph"),
  ctx = graph.getContext("2d");
function status(text, error = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", error);
}
function options() {
  const opts = Object.fromEntries(
    ["k", "support", "mapq", "baseq"].map((id) => [id, +$(id).value]),
  );
  opts.minVaf = +$("minVaf").value / 100;
  if (!Number.isFinite(opts.minVaf) || opts.minVaf < 0 || opts.minVaf > 1)
    throw Error("Invalid VAF threshold");
  opts.duplicates = $("duplicates").checked;
  if (
    !Number.isInteger(opts.support) ||
    opts.support < 1 ||
    opts.support > 1000 ||
    !Number.isInteger(opts.mapq) ||
    opts.mapq < 0 ||
    opts.mapq > 254 ||
    !Number.isInteger(opts.baseq) ||
    opts.baseq < 0 ||
    opts.baseq > 93
  )
    throw Error("Check evidence thresholds");
  return opts;
}
function updateButtons() {
  for (const id of ["open", "demo", "find"]) $(id).disabled = busy;
  $("run").disabled = busy || !opened;
  for (const id of ["query", "gfa", "allreads"])
    $(id).disabled = busy || !result;
  $("expand").disabled = busy || !selection;
}
function post(m) {
  if (busy) return;
  busy = true;
  updateButtons();
  timer = setTimeout(() => {
    destroyWorker();
    opened = false;
    busy = false;
    status(
      "Job exceeded 90 seconds and was stopped. Narrow the interval, then reopen the BAM. No reads were downsampled.",
      true,
    );
    updateButtons();
  }, 90000);
  worker.postMessage(m);
}
function destroyWorker() {
  clearTimeout(timer);
  worker?.terminate();
  if (workerURL) URL.revokeObjectURL(workerURL);
  worker = null;
  workerURL = null;
}
function startWorker() {
  destroyWorker();
  workerURL = URL.createObjectURL(
    new Blob([WORKER_SOURCE], { type: "text/javascript" }),
  );
  worker = new Worker(workerURL);
  worker.onmessage = receive;
  worker.onerror = (e) => {
    clearTimeout(timer);
    busy = false;
    status("Worker error: " + e.message, true);
    updateButtons();
  };
}
function receive({ data: m }) {
  if (m.type === "progress") {
    status(m.message);
    return;
  }
  clearTimeout(timer);
  busy = false;
  if (m.type === "error") status(m.message, true);
  if (m.type === "opened") {
    opened = true;
    contigs = m.contigs;
    $("contigs").replaceChildren(
      ...contigs.map((c) => {
        const o = make("option");
        o.value = c.refName;
        return o;
      }),
    );
    status(
      `BAM opened · ${contigs.length} contigs · ${m.reference ? "reference annotation available" : "reference-free assembly"}. Choose a local interval.`,
    );
  }
  if (m.type === "result") {
    const savedReads =
      m.preserveView && readsOpen
        ? readRows.flatMap((rr) => rr.map((r) => r.id))
        : null;
    const savedReadView = readView,
      savedScroll = $("reads").scrollTop;
    result = m.result;
    selection = null;
    query = null;
    queried.clear();
    highlight.clear();
    $("linkage").hidden = true;
    $("evidence").hidden = true;
    readsOpen = false;
    readView = null;
    readRows = [];
    readBaseline.clear();
    readCanvas = null;
    readCtx = null;
    $("reads").replaceChildren();
    $("chr").value = result.region.chr;
    $("start").value = result.region.start + 1;
    $("end").value = result.region.end;
    const s = result.stats;
    status(
      `${s.records} input records retained · ${s.eligibleRecords} quality-eligible records · ${s.fragments} physical fragments assembled. No downsampling.`,
    );
    $("metrics").replaceChildren(
      ...[
        [s.fragments, "eligible fragments"],
        [result.unitigs.length, "sequence segments"],
        [result.bubbles.length, "observed bubbles"],
        [result.blocks.length, "local evidence blocks"],
      ].map(([n, t]) => {
        const e = make("span");
        e.append(make("strong", n), document.createTextNode(t));
        return e;
      }),
    );
    $("graph-note").textContent =
      `${result.region.chr}:${result.region.start + 1}–${result.region.end} · topology layout, not a genomic axis · ${s.orientedKmers.toLocaleString()} oriented k-mers · ${s.weakUnitigs} weak segments retained · ${s.complex} bounded/complex fork searches unresolved. ${
        Object.entries(s.issues)
          .map(([k, v]) => `${v} ${k}`)
          .join("; ") || "No records excluded by quality/geometry filters."
      }`;
    layout();
    renderBlocks();
    renderStructural();
    renderSegments();
    inspect(null);
    if (m.preserveView) draw();
    else fit();
    if (savedReads) {
      readView = savedReadView;
      showReadIDs(savedReads);
      $("reads").scrollTop = savedScroll;
      drawReads();
    }
  }
  if (m.type === "query") {
    query = m.query;
    renderQuery();
    markQueries();
    status(
      "Variant linkage evaluated from physical fragments. Gold graph outlines mark matching sequence evidence.",
    );
  }
  if (m.type === "gfa")
    download(m.text, "phase-local-assembly.gfa", "text/plain");
  updateButtons();
}
function download(text, name, type) {
  const url = URL.createObjectURL(new Blob([text], { type })),
    a = make("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
$("open").onclick = () => {
  try {
    const files = Object.fromEntries(
      ["bam", "index", "reference", "fai"].map((id) => [id, $(id).files[0]]),
    );
    if (!files.bam || !files.index)
      throw Error("Choose a BAM and matching BAI / CSI");
    startWorker();
    opened = false;
    clearResult();
    status("Opening local BAM header and index…");
    post({ type: "open", files });
  } catch (e) {
    status(e.message, true);
  }
};
watchVaf($("minVaf"), {
  current: () => result?.settings.minVaf,
  busy: () => busy,
  update: (minVaf) => {
    status("Updating visible variants from loaded reads…");
    post({ type: "vaf", minVaf });
  },
  error: (e) => status(e.message, true),
});
$("demo").onclick = () => {
  try {
    if (!worker) startWorker();
    status("Assembling synthetic paired reads…");
    post({
      type: "demo",
      scenario: $("demo-kind").value,
      opts: options(),
      reference: $("demo-ref").checked,
    });
  } catch (e) {
    status(e.message, true);
  }
};
$("run").onclick = () => {
  try {
    const region = {
      chr: $("chr").value.trim(),
      start: +$("start").value - 1,
      end: +$("end").value,
    };
    status("Preparing local sequence assembly…");
    post({ type: "analyze", region, opts: options() });
  } catch (e) {
    status(e.message, true);
  }
};
$("gfa").onclick = () => post({ type: "gfa" });
function loadCatalog() {
  const b = $("build").value;
  if (b !== catalogBuild) {
    catalog = JSON.parse(
      new TextDecoder().decode(
        ungzip(
          Uint8Array.from(atob(b === "hg19" ? hg19 : hg38), (c) =>
            c.charCodeAt(0),
          ),
        ),
      ),
    );
    catalogBuild = b;
  }
  return catalog;
}
$("find").onclick = () => {
  try {
    const text = $("gene").value.trim().toUpperCase();
    if (!text) throw Error("Enter a gene symbol");
    genes = loadCatalog().filter((g) => g[0].toUpperCase() === text);
    if (!genes.length)
      throw Error("Gene not found in the selected offline catalog");
    $("genes").replaceChildren(
      ...genes.map((g, i) => {
        const o = make("option", `${g[0]} · ${g[1]}:${g[3] + 1}–${g[4]}`);
        o.value = i;
        return o;
      }),
    );
    chooseGene();
  } catch (e) {
    status(e.message, true);
  }
};
function chooseGene() {
  const g = genes[+$("genes").value];
  if (!g) return;
  $("transcript").replaceChildren(
    ...g[6].map((t, i) => {
      const o = make("option", t[0]);
      o.value = i;
      return o;
    }),
  );
  chooseTranscript();
}
function chooseTranscript() {
  const g = genes[+$("genes").value],
    t = g?.[6][+$("transcript").value];
  if (!t) return;
  selected = { chr: g[1], strand: g[2], transcript: t };
  const aliases = [g[1], g[1].startsWith("chr") ? g[1].slice(3) : "chr" + g[1]];
  $("chr").value =
    contigs.find((c) => aliases.includes(c.refName))?.refName || g[1];
  $("start").value = t[1] + 1;
  $("end").value = t[2];
  $("exon").replaceChildren(
    make("option", "Whole transcript interval"),
    ...t[3].map(([s, e], i) => {
      const o = make(
        "option",
        `Exon ${selected.strand === "-" ? t[3].length - i : i + 1} · ${s + 1}–${e}`,
      );
      o.value = i;
      return o;
    }),
  );
  $("exon").firstChild.value = "";
  status(
    `Selected ${g[0]} · ${t[0]} · ${$("build").value}. Check that your BAM uses this assembly, then assemble.`,
  );
}
$("genes").onchange = chooseGene;
$("transcript").onchange = chooseTranscript;
$("exon").onchange = () => {
  const g = genes[+$("genes").value],
    t = g?.[6][+$("transcript").value];
  if (!t) return;
  const interval =
    $("exon").value === "" ? [t[1], t[2]] : t[3][+$("exon").value];
  $("start").value = interval[0] + 1;
  $("end").value = interval[1];
  status(
    "Selected target interval updated. Assemble to rebuild this local region.",
  );
};
$("build").onchange = () => {
  catalog = null;
  selected = null;
  genes = [];
  $("genes").replaceChildren(make("option", "Choose a gene"));
  $("transcript").replaceChildren(make("option", "Choose a transcript"));
  const o = make("option", "Whole transcript interval");
  o.value = "";
  $("exon").replaceChildren(o);
};
$("query").onclick = () => {
  try {
    post({
      type: "query",
      text: $("variants").value,
      catalog: loadCatalog(),
      selected,
    });
  } catch (e) {
    status(e.message, true);
  }
};
function clearResult() {
  result = null;
  selection = null;
  query = null;
  readFragments = [];
  readView = null;
  highlight.clear();
  queried.clear();
  positions = [];
  readRows = [];
  readBaseline.clear();
  readCanvas = null;
  readCtx = null;
  readDrag = null;
  readsOpen = false;
  $("blocks").replaceChildren(
    make(
      "div",
      "Assembly branches and their linkage will appear here.",
      "empty",
    ),
  );
  $("segments").replaceChildren();
  $("structural").replaceChildren();
  $("structural-summary").textContent = "Structural evidence";
  $("graph-note").textContent =
    "Topology layout, not a genomic coordinate axis. Both strands are folded into oriented sequence segments.";
  $("metrics").replaceChildren(
    make("span", "Reference optional"),
    make("span", "No fixed ploidy"),
    make("span", "No read downsampling"),
  );
  $("evidence").hidden = true;
  $("linkage").hidden = true;
  $("reads").replaceChildren();
  inspect(null);
  draw();
}
$("reset").onclick = () => {
  destroyWorker();
  busy = false;
  opened = false;
  contigs = [];
  catalog = null;
  catalogBuild = null;
  genes = [];
  selected = null;
  for (const id of [
    "bam",
    "index",
    "reference",
    "fai",
    "gene",
    "chr",
    "start",
    "end",
    "variants",
  ])
    $(id).value = "";
  $("contigs").replaceChildren();
  $("genes").replaceChildren(make("option", "Choose a gene"));
  $("transcript").replaceChildren(make("option", "Choose a transcript"));
  const o = make("option", "Whole transcript interval");
  o.value = "";
  $("exon").replaceChildren(o);
  clearResult();
  startWorker();
  status("Session cleared. Case data released; no files were saved.");
  updateButtons();
};

function layout() {
  if (!result) return;
  positions = new Array(result.unitigs.length);
  const adjacent = result.unitigs.map(() => new Set());
  for (const e of result.links) {
    adjacent[e.from].add(e.to);
    adjacent[e.to].add(e.from);
  }
  const seen = new Set();
  let yOffset = 0;
  for (const start of result.unitigs
    .slice()
    .sort((a, b) => a.position - b.position)) {
    if (seen.has(start.id)) continue;
    const component = [start.id];
    seen.add(start.id);
    for (let i = 0; i < component.length; i++)
      for (const id of adjacent[component[i]])
        if (!seen.has(id)) {
          seen.add(id);
          component.push(id);
        }
    const sorted = component
        .map((id) => result.unitigs[id])
        .sort((a, b) => a.position - b.position || a.id - b.id),
      lanes = [];
    let column = 0,
      anchor = sorted[0].position;
    for (const u of sorted) {
      if (u.position - anchor > result.settings.k) {
        column++;
        anchor = u.position;
      }
      const x = column * 205,
        w = Math.min(165, 80 + Math.log2(u.length) * 9);
      let lane = lanes.findIndex((end) => end < x - 20);
      if (lane < 0) lane = lanes.length;
      lanes[lane] = x + w;
      positions[u.id] = { x, y: yOffset + lane * 65, w, h: 30 };
    }
    yOffset += lanes.length * 65 + 65;
  }
}
function resizeCanvas(canvas, height) {
  const ratio = devicePixelRatio || 1,
    w = canvas.clientWidth || 800;
  canvas.width = Math.round(w * ratio);
  canvas.height = Math.round(height * ratio);
  const c = canvas.getContext("2d");
  c.setTransform(ratio, 0, 0, ratio, 0, 0);
  return { c, w, h: height };
}
function fit() {
  if (!positions.length) {
    draw();
    return;
  }
  const maxX = positions.reduce((max, p) => Math.max(max, p.x + p.w), 0),
    maxY = positions.reduce((max, p) => Math.max(max, p.y + p.h), 0),
    w = graph.clientWidth || 800;
  camera.scale = Math.min(1.4, (w - 50) / (maxX + 20), 380 / (maxY + 30));
  camera.x = (w - maxX * camera.scale) / 2;
  camera.y = (440 - maxY * camera.scale) / 2;
  draw();
}
function draw() {
  if (result) graph.dataset.minVaf = result.settings.minVaf;
  const { c, w, h } = resizeCanvas(graph, 440);
  c.fillStyle = "#fbfdfd";
  c.fillRect(0, 0, w, h);
  c.strokeStyle = "#edf2f3";
  c.lineWidth = 1;
  for (let x = 12; x < w; x += 24)
    for (let y = 12; y < h; y += 24) {
      c.beginPath();
      c.moveTo(x, y);
      c.lineTo(x + 1, y);
      c.stroke();
    }
  if (!result) {
    c.fillStyle = "#8b9ba1";
    c.font = "13px -apple-system,sans-serif";
    c.textAlign = "center";
    c.fillText("Assemble a locus to reveal its sequence graph", w / 2, h / 2);
    return;
  }
  c.save();
  c.translate(camera.x, camera.y);
  c.scale(camera.scale, camera.scale);
  for (const e of result.links) {
    const a = positions[e.from],
      b = positions[e.to];
    if (!a || !b) continue;
    const x = a.x + (e.fromOrientation === "+" ? a.w : 0),
      y = a.y + 15,
      x2 = b.x + (e.toOrientation === "+" ? 0 : b.w),
      y2 = b.y + 15,
      active = highlight.has(e.from) && highlight.has(e.to);
    c.strokeStyle = active ? "#df9b37" : "#b8cbd1";
    c.lineWidth = active ? 3 : 1.5;
    c.setLineDash(e.support < result.settings.support ? [5, 4] : []);
    c.beginPath();
    c.moveTo(x, y);
    const out = e.fromOrientation === "+" ? 1 : -1,
      into = e.toOrientation === "+" ? 1 : -1;
    if (e.from === e.to)
      c.bezierCurveTo(x + out * 65, y - 55, x2 - into * 65, y2 - 55, x2, y2);
    else c.bezierCurveTo(x + out * 55, y, x2 - into * 55, y2, x2, y2);
    c.stroke();
    c.setLineDash([]);
    c.fillStyle = c.strokeStyle;
    c.beginPath();
    c.moveTo(x2, y2);
    c.lineTo(x2 - into * 6, y2 - 3);
    c.lineTo(x2 - into * 6, y2 + 3);
    c.fill();
  }
  for (const u of result.unitigs) {
    const p = positions[u.id];
    const active = highlight.has(u.id),
      gold = queried.has(u.id);
    c.globalAlpha = highlight.size && !active ? 0.35 : 1;
    c.fillStyle =
      u.reference !== null && u.reference > 0.99 ? "#deefeb" : "#eee7f6";
    c.strokeStyle = gold
      ? "#d89325"
      : active
        ? "#15887d"
        : u.support < result.settings.support
          ? "#b5c4c9"
          : u.reference !== null && u.reference > 0.99
            ? "#398f85"
            : "#9878bb";
    c.lineWidth = gold || active ? 3 : 1.4;
    c.setLineDash(u.support < result.settings.support ? [4, 3] : []);
    c.beginPath();
    c.roundRect(p.x, p.y, p.w, p.h, 7);
    c.fill();
    c.stroke();
    c.setLineDash([]);
    if (p.w * camera.scale < 36) {
      c.globalAlpha = 1;
      continue;
    }
    c.fillStyle = "#344f58";
    c.textAlign = "left";
    c.font = `${Math.max(11, 8 / camera.scale)}px -apple-system,sans-serif`;
    c.fillText(
      `u${u.id + 1} · ${u.length}${camera.scale < 0.7 ? "" : " bp"}`,
      p.x + 9,
      p.y + 19,
    );
    c.fillStyle = "#7c8f96";
    c.font = `${Math.max(9, 7 / camera.scale)}px -apple-system,sans-serif`;
    c.fillText(
      camera.scale < 0.7
        ? `${u.support} frags${u.cyclic ? " · loop" : ""}`
        : `${u.support} min. fragment support${u.cyclic ? " · loop" : ""}`,
      p.x + 2,
      p.y + 44,
    );
    c.globalAlpha = 1;
  }
  c.restore();
}
$("fit").onclick = fit;
function zoom(f, x = graph.clientWidth / 2, y = 220) {
  const before = camera.scale,
    next = Math.max(0.02, Math.min(6, before * f));
  camera.x = x - ((x - camera.x) * next) / before;
  camera.y = y - ((y - camera.y) * next) / before;
  camera.scale = next;
  draw();
}
$("zoomin").onclick = () => zoom(1.25);
$("zoomout").onclick = () => zoom(0.8);
graph.onwheel = (e) => {
  e.preventDefault();
  const r = graph.getBoundingClientRect();
  zoom(Math.exp(-e.deltaY * 0.001), e.clientX - r.left, e.clientY - r.top);
};
let drag = null;
graph.onpointerdown = (e) => {
  graph.setPointerCapture(e.pointerId);
  drag = {
    x: e.clientX,
    y: e.clientY,
    cx: camera.x,
    cy: camera.y,
    moved: false,
  };
};
graph.onpointermove = (e) => {
  if (!drag) return;
  const dx = e.clientX - drag.x,
    dy = e.clientY - drag.y;
  if (Math.abs(dx) + Math.abs(dy) > 4) drag.moved = true;
  camera.x = drag.cx + dx;
  camera.y = drag.cy + dy;
  draw();
};
graph.onpointerup = (e) => {
  if (!drag) return;
  if (!drag.moved && result) {
    const r = graph.getBoundingClientRect(),
      x = (e.clientX - r.left - camera.x) / camera.scale,
      y = (e.clientY - r.top - camera.y) / camera.scale;
    const id = positions.findIndex(
      (p) => x >= p.x && x <= p.x + p.w && y >= p.y && y <= p.y + p.h,
    );
    if (id >= 0) selectUnitig(id);
  }
  drag = null;
};
graph.onpointercancel = () => (drag = null);
graph.onkeydown = (e) => {
  if (["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"].includes(e.key)) {
    e.preventDefault();
    camera.x += e.key === "ArrowLeft" ? 40 : e.key === "ArrowRight" ? -40 : 0;
    camera.y += e.key === "ArrowUp" ? 40 : e.key === "ArrowDown" ? -40 : 0;
    draw();
  }
};
new ResizeObserver(() => {
  draw();
  if (readsOpen) drawReads();
}).observe(graph);

function renderSegments(page = 0) {
  const root = $("segments");
  root.replaceChildren();
  const start = page * 100,
    end = Math.min(start + 100, result.unitigs.length);
  if (result.unitigs.length > 100) {
    root.append(
      make("span", `${start + 1}–${end} of ${result.unitigs.length}`),
    );
    for (const [text, step] of [
      ["Previous 100", -1],
      ["Next 100", 1],
    ]) {
      const button = make("button", text);
      button.disabled = step < 0 ? page === 0 : end === result.unitigs.length;
      button.onclick = () => renderSegments(page + step);
      root.append(button);
    }
  }
  for (const u of result.unitigs.slice(start, end)) {
    const b = make(
      "button",
      `u${u.id + 1}: ${u.length} bp · ${u.support} support`,
    );
    b.onclick = () => {
      const p = positions[u.id];
      camera.scale = Math.max(0.9, camera.scale);
      camera.x = graph.clientWidth / 2 - (p.x + p.w / 2) * camera.scale;
      camera.y = 220 - (p.y + 15) * camera.scale;
      selectUnitig(u.id);
    };
    root.append(b);
  }
}
function renderBlocks() {
  const root = $("blocks");
  root.replaceChildren();
  if (!result.blocks.length) {
    root.append(
      make(
        "p",
        "No multi-allelic, read-spanned bubbles resolved. A single graph path does not establish one cellular clone. Inspect the graph and reads for coverage gaps, weak branches and repeat ambiguity.",
        "note",
      ),
    );
  }
  for (const b of result.blocks) {
    const box = make("div", undefined, "block"),
      n = b.groups.filter((g) => g.supported).length;
    box.append(
      make(
        "div",
        `Block ${b.id + 1} · ${n} supported local pattern${n === 1 ? "" : "s"}`,
        "block-title",
      ),
      make(
        "p",
        `${b.total} informative fragments · ${b.features.length} linked bubbles · ${b.shared.length} shared/ambiguous · ${b.contradictory.length} internally contradictory. Weak patterns stay separate.`,
        "block-note",
      ),
    );
    const weak = make("details"),
      weakBody = make("div");
    const weakCount = b.groups.filter((g) => !g.supported).length;
    weak.append(
      make("summary", `${weakCount} weak patterns · expand to inspect`),
      weakBody,
    );
    for (const g of b.groups) {
      const button = make(
        "button",
        undefined,
        "path" + (g.supported ? "" : " weak"),
      );
      button.dataset.block = b.id;
      button.dataset.group = g.id;
      button.append(
        make("strong", `H${b.id + 1}.${g.id + 1}`),
        make(
          "span",
          `${g.fragments.length}/${b.total} (${((100 * g.fragments.length) / b.total).toFixed(1)}%)${g.supported ? "" : " · weak"}`,
        ),
      );
      const pattern = make("span", undefined, "pattern");
      for (const id of b.features) {
        const a = new Map(g.calls).get(id),
          chip = make(
            "span",
            a === undefined ? "?" : `${id + 1}:${String.fromCharCode(65 + a)}`,
            "allele" + (a === undefined ? "" : ` alt-${a % 4}`),
          );
        chip.title = `Bubble ${id + 1}, allele ${a === undefined ? "unknown" : String.fromCharCode(65 + a)}`;
        pattern.append(chip);
      }
      button.append(pattern);
      button.onclick = () => selectGroup(b.id, g.id);
      (g.supported ? box : weakBody).append(button);
    }
    if (weakCount) box.append(weak);
    root.append(box);
  }
  root.append(
    make(
      "p",
      `${result.stats.unthreaded} eligible fragments have no resolved bubble allele; they remain in the graph/read evidence. Fractions above are exclusive assignments among informative fragments, not VAF or cellular abundance.`,
      "block-note",
    ),
  );
}
function inspect(s) {
  selection = s;
  $("selection-title").textContent = s?.title || "Explore the reconstruction";
  $("inspection").replaceChildren();
  $("expand").disabled = !s;
  if (!s) {
    $("inspection").append(
      make(
        "p",
        "Select a graph segment or suggested haplotype to inspect sequence and paired reads.",
        "muted",
      ),
    );
    return;
  }
  s.render($("inspection"));
}
function selectUnitig(id) {
  for (const el of document.querySelectorAll(".path.selected"))
    el.classList.remove("selected");
  const u = result.unitigs[id];
  highlight = new Set([id]);
  inspect({
    title: `Sequence segment u${id + 1}`,
    fragments: u.fragments,
    render: (root) => {
      root.append(
        make(
          "p",
          `${u.length} bp · minimum support ${u.support} distinct fragments per k-mer`,
          "fact",
        ),
        make(
          "p",
          `${u.fragments.length} fragments traverse at least part of this segment. ${u.reference === null ? "Reference annotation unavailable." : `${Math.round(u.reference * 100)}% of constituent k-mers occur in the reference interval.`}`,
          "muted",
        ),
      );
      if (u.cyclic)
        root.append(
          make(
            "p",
            "This segment participates in a cycle. Repeat multiplicity and copy number are unresolved.",
            "note",
          ),
        );
      root.append(make("div", u.seq, "sequence"));
    },
  });
  draw();
  if (readsOpen) showReads(u.fragments);
}
function selectGroup(blockID, groupID) {
  const b = result.blocks[blockID],
    g = b.groups[groupID],
    fs = new Set(g.fragments);
  highlight = new Set(
    result.unitigs
      .filter((u) => u.fragments.some((id) => fs.has(id)))
      .map((u) => u.id),
  );
  for (const el of document.querySelectorAll(".path"))
    el.classList.toggle(
      "selected",
      +el.dataset.block === blockID && +el.dataset.group === groupID,
    );
  inspect({
    title: `Local pattern H${blockID + 1}.${groupID + 1}`,
    fragments: g.fragments,
    render: (root) => {
      root.append(
        make(
          "p",
          `${g.fragments.length} exclusively assigned fragments / ${b.total} informative (${((100 * g.fragments.length) / b.total).toFixed(1)}%)`,
          "fact",
        ),
        make(
          "p",
          `${b.shared.filter((s) => s.groups.includes(g.id)).length} additional fragments compatible with this and other patterns are shared, not counted twice.`,
          "muted",
        ),
      );
      for (const [id, a] of g.calls) {
        const bubble = result.bubbles[id],
          allele = bubble.alleles[a];
        root.append(
          make(
            "p",
            `Bubble ${id + 1} · allele ${String.fromCharCode(65 + a)} · ${allele.callSupport}/${bubble.callable} callable fragments (${bubble.callable ? ((100 * allele.callSupport) / bubble.callable).toFixed(1) : "0"}%)${allele.isReference === true ? " · reference" : allele.isReference === false ? " · differs from reference" : ""}`,
            "fact",
          ),
        );
        const seq = make("div", undefined, "sequence-diff");
        if (allele.isReference) {
          seq.textContent =
            "Matches the supplied reference across this bubble.";
          root.append(seq);
        } else if (bubble.referenceSequence) {
          const ref = bubble.referenceSequence,
            alt = allele.seq;
          let p = 0,
            z = 0;
          while (p < ref.length && p < alt.length && ref[p] === alt[p]) p++;
          while (
            z < ref.length - p &&
            z < alt.length - p &&
            ref.at(-1 - z) === alt.at(-1 - z)
          )
            z++;
          seq.append(
            document.createTextNode(alt.slice(Math.max(0, p - 10), p)),
            make("mark", alt.slice(p, alt.length - z) || "∅"),
            document.createTextNode(
              alt.slice(
                alt.length - z,
                Math.min(alt.length, alt.length - z + 10),
              ),
            ),
          );
          root.append(
            seq,
            make(
              "p",
              `Anchor-trimmed reference: ${ref.slice(p, ref.length - z) || "∅"} · allele: ${alt.slice(p, alt.length - z) || "∅"}`,
              "mini-legend",
            ),
          );
        } else {
          seq.textContent = allele.seq;
          root.append(seq);
        }
      }
      const distances = b.groups
        .filter((h) => h.id !== g.id)
        .map((h) => {
          const other = new Map(h.calls);
          let same = 0,
            diff = 0;
          for (const [id, a] of g.calls)
            if (other.has(id)) {
              if (other.get(id) === a) same++;
              else diff++;
            }
          return { h, same, diff };
        })
        .filter((d) => d.same + d.diff)
        .sort(
          (a, b) => a.diff / (a.same + a.diff) - b.diff / (b.same + b.diff),
        );
      if (distances.length) {
        const d = distances[0];
        root.append(
          make(
            "p",
            `Closest observed allele pattern: H${blockID + 1}.${d.h.id + 1}; ${d.same} shared and ${d.diff} differing bubble alleles. Similarity does not establish ancestry. Missing alleles are ignored.`,
            "note",
          ),
        );
      }
      root.append(
        make(
          "p",
          "Compatible partial observations suggest this local pattern. Long-range phase, repeat copy number and cellular clone abundance remain unresolved. No fixed number of haplotypes is assumed.",
          "note",
        ),
      );
    },
  });
  draw();
  if (readsOpen) showReads(g.fragments);
}
function renderStructural() {
  const s = result.structural;
  $("structural-summary").textContent =
    `Structural evidence · ${s.events.length} local split/mate clues`;
  const root = $("structural");
  root.replaceChildren(
    make(
      "p",
      "These are local alignment clues, not confirmed rearrangements. Loops/branches may reflect repeats, duplications, deletions or errors; a graph alone does not uniquely identify inversion or translocation partners.",
      "muted",
    ),
  );
  for (const e of s.events) {
    const item = make("div", undefined, "event");
    item.append(
      make("strong", e.type),
      make(
        "p",
        `${e.local.chr}:${e.local.pos + 1} (${e.local.strand}) → ${e.partner.chr}:${e.partner.pos + 1} (${e.partner.strand}) · ${e.fragments} fragments; ${e.eligibleFragments} pass both-end quality checks.`,
      ),
    );
    for (const n of e.notes) item.append(make("p", n, "muted"));
    const b = make("button", "Inspect supporting reads");
    b.onclick = () => showReadIDs(e.records);
    item.append(b);
    root.append(item);
  }
  for (const w of s.warnings) root.append(make("p", w, "note"));
}

function renderQuery() {
  const root = $("query-results");
  root.replaceChildren();
  $("linkage").hidden = false;
  for (const e of query.entries)
    root.append(
      make(
        "div",
        `${e.line}. ${e.input} — ${e.error || `${result.region.chr}:${e.pos + 1} ${e.ref || "∅"}>${e.alt || "∅"}${e.unverified ? " · REF unverified" : ""}`}`,
        "query-item" + (e.error ? " error" : ""),
      ),
    );
  for (const p of query.pairs) {
    const item = make("div", undefined, "pair-result");
    item.append(
      make("strong", p.status),
      make("div", p.left ? `${p.left} ↔ ${p.right}` : `Lines ${p.a} ↔ ${p.b}`),
      make(
        "div",
        `Both ALT: ${p.both} · left ALT/right REF: ${p.leftOnly} · left REF/right ALT: ${p.rightOnly} · both REF: ${p.neither}`,
      ),
      make("div", `Unknown: ${p.unknown} · contradictory/other: ${p.other}`),
    );
    root.append(item);
  }
  root.append(
    make(
      "p",
      "A linked pair may span an unsequenced insert; linkage assumes the BAM query name/read group identifies one physical fragment. Separate patterns do not establish trans or a fixed ploidy.",
      "note",
    ),
  );
  if (readsOpen) drawReads();
}
function markQueries() {
  queried.clear();
  for (const entry of query.entries) {
    if (!entry.variant) continue;
    for (const r of result.reads) {
      if (r.exclusion || !hasSequence(r)) continue;
      for (const run of alignmentRuns(r)) {
        if (
          !"M=X".includes(run.op) ||
          entry.pos < run.refStart ||
          entry.pos >= run.refStart + run.length
        )
          continue;
        const q = run.queryStart + entry.pos - run.refStart,
          half = Math.floor(result.settings.k / 2),
          key = r.seq.slice(q - half, q - half + result.settings.k);
        if (key.length !== result.settings.k) continue;
        const rev = reverseComplement(key);
        for (const u of result.unitigs)
          if (u.seq.includes(key) || u.seq.includes(rev)) queried.add(u.id);
      }
    }
  }
  draw();
}

// Read rows are virtualized for rendering only. All fragment rows remain
// reachable by the native vertical scrollbar; evidence is never sampled.
let readCanvas,
  readCtx,
  readRows = [],
  readBaseline = new Map(),
  readDrag = null;
function showReadIDs(ids) {
  const groups = new Map();
  for (const id of ids) {
    const r = result.reads[id];
    if (!r) continue;
    const key = JSON.stringify([r.sample, r.rg, r.name]);
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(r);
  }
  openRows([...groups.values()]);
}
function showReads(ids) {
  readFragments = ids;
  const all = new Set(ids.flatMap((id) => result.fragments[id]?.reads || []));
  showReadIDs([...all]);
}
function openRows(rows) {
  readRows = rows;
  readsOpen = true;
  $("evidence").hidden = false;
  $("read-count").textContent =
    `${rows.length} fragments · ${rows.reduce((n, r) => n + r.length, 0)} alignment records · all rows available by vertical scrolling`;
  if (!readView)
    readView = { start: result.region.start, end: result.region.end };
  readBaseline = new Map();
  const counts = new Map();
  for (const f of result.fragments)
    for (const [p, b] of f.calls) {
      if (!counts.has(p)) counts.set(p, new Map());
      const c = counts.get(p);
      c.set(b, (c.get(b) || 0) + 1);
    }
  for (const [p, c] of counts)
    readBaseline.set(p, [...c].sort((a, b) => b[1] - a[1])[0][0]);
  const space = make("div");
  space.id = "read-space";
  readCanvas = make("canvas");
  readCanvas.id = "read-canvas";
  readCanvas.height = 360;
  readCtx = readCanvas.getContext("2d");
  space.append(readCanvas);
  $("reads").replaceChildren(space);
  const rule = [...document.styleSheets[0].cssRules].find(
    (r) => r.selectorText === "#read-space",
  );
  rule.style.height = Math.max(360, rows.length * 49 + 30) + "px";
  $("reads").scrollTop = 0;
  $("reads").onscroll = drawReads;
  readCanvas.onpointerdown = (e) => {
    readCanvas.setPointerCapture(e.pointerId);
    readDrag = { x: e.clientX, start: readView.start, end: readView.end };
  };
  readCanvas.onpointermove = (e) => {
    if (!readDrag) return;
    const w = readCanvas.clientWidth - 180,
      shift = ((e.clientX - readDrag.x) / w) * (readDrag.end - readDrag.start);
    const span = readDrag.end - readDrag.start;
    readView.start = Math.max(0, readDrag.start - shift);
    readView.end = readView.start + span;
    drawReads();
  };
  readCanvas.onpointerup = readCanvas.onpointercancel = () => (readDrag = null);
  readCanvas.onwheel = (e) => {
    if (!e.ctrlKey && !e.metaKey) return;
    e.preventDefault();
    const center = (readView.start + readView.end) / 2,
      span = Math.max(
        40,
        Math.min(
          1000000,
          (readView.end - readView.start) * Math.exp(e.deltaY * 0.002),
        ),
      );
    readView = {
      start: Math.max(0, center - span / 2),
      end: Math.max(0, center - span / 2) + span,
    };
    drawReads();
  };
  drawReads();
}
function drawReads() {
  if (!readCanvas || !result) return;
  const { c, w } = resizeCanvas(readCanvas, 360),
    scroll = $("reads").scrollTop,
    first = Math.max(0, Math.floor((scroll - 30) / 49)),
    last = Math.min(readRows.length, first + 10),
    left = 180,
    span = readView.end - readView.start,
    x = (p) => left + ((p - readView.start) / span) * (w - left - 12);
  c.fillStyle = "#fbfdfd";
  c.fillRect(0, 0, w, 360);
  c.font = "10px -apple-system,sans-serif";
  c.fillStyle = "#7b9098";
  c.fillText(
    `${result.region.chr}:${Math.floor(readView.start) + 1}–${Math.ceil(readView.end)} · ctrl+wheel to zoom`,
    left,
    14,
  );
  c.save();
  c.beginPath();
  c.rect(left, 20, w - left, 340);
  c.clip();
  for (const e of query?.entries || [])
    if (e.variant) {
      c.fillStyle = "#e8b95630";
      c.fillRect(x(e.pos) - 2, 20, 4, 340);
    }
  for (let i = first; i < last; i++) {
    const rs = readRows[i],
      y = 30 + i * 49 - scroll + 20;
    const starts = rs.filter((r) => r.alignmentUsable).map((r) => r.start),
      ends = rs.filter((r) => r.alignmentUsable).map((r) => r.end);
    if (starts.length > 1) {
      c.strokeStyle = "#bdcbd0";
      c.setLineDash([3, 3]);
      c.beginPath();
      c.moveTo(x(Math.min(...starts)), y);
      c.lineTo(x(Math.max(...ends)), y);
      c.stroke();
      c.setLineDash([]);
    }
    for (const r of rs) {
      const yy = y + (r.flags & 128 ? 7 : -3);
      for (const run of alignmentRuns(r)) {
        const { refStart: p, queryStart: q, length: n, op } = run;
        if ("M=X".includes(op)) {
          c.fillStyle = r.exclusion
            ? "#cbd2d4"
            : r.flags & 128
              ? "#a9c7db"
              : "#93beb6";
          if (hasSequence(r))
            c.fillRect(x(p), yy, Math.max(1, x(p + n) - x(p)), 6);
          else {
            c.strokeStyle = "#aab9bd";
            c.strokeRect(x(p), yy, Math.max(1, x(p + n) - x(p)), 6);
          }
          if (hasSequence(r))
            for (let j = 0; j < n; j++) {
              const pos = p + j,
                base = r.seq[q + j],
                ref =
                  result.reference?.[pos - result.region.start] ||
                  readBaseline.get(pos);
              if (base && ref && base !== ref && !hiddenBase(r, q + j)) {
                c.fillStyle =
                  { A: "#289964", C: "#3987ce", G: "#d79526", T: "#d25a69" }[
                    base
                  ] || "#777";
                c.fillRect(x(pos), yy, Math.max(2, x(pos + 1) - x(pos)), 7);
                if (span < 180) {
                  c.fillStyle = "#344f58";
                  c.fillText(base, x(pos), yy - 2);
                }
              }
            }
        } else if (op === "I" && !hiddenEvent(r, op, p)) {
          c.fillStyle = "#9266bd";
          c.beginPath();
          c.moveTo(x(p), yy - 4);
          c.lineTo(x(p) - 4, yy - 10);
          c.lineTo(x(p) + 4, yy - 10);
          c.fill();
        } else if ((op === "D" || op === "N") && !hiddenEvent(r, op, p)) {
          c.strokeStyle = "#c66d76";
          c.beginPath();
          c.moveTo(x(p), yy + 3);
          c.lineTo(x(p + n), yy + 3);
          c.stroke();
        } else if (op === "S") {
          c.fillStyle = "#c2b2d3";
          c.fillRect(x(p) - 2, yy, 3, 7);
        }
      }
    }
  }
  c.restore();
  c.fillStyle = "#fbfdfd";
  c.fillRect(0, 0, left - 4, 360);
  for (let i = first; i < last; i++) {
    const rs = readRows[i],
      y = 30 + i * 49 - scroll + 20;
    c.fillStyle = "#486774";
    c.fillText(rs[0].name.slice(0, 26), 7, y + 3);
    c.fillStyle = "#8b9da3";
    c.font = "9px -apple-system,sans-serif";
    c.fillText(
      rs
        .map(
          (r) =>
            `${r.flags & 128 ? "R2" : r.flags & 64 ? "R1" : "read"} ${r.displayIssue || r.cigar || "no CIGAR"}${r.exclusion ? " [excluded]" : ""}`,
        )
        .join(" · ")
        .slice(0, 40),
      7,
      y + 17,
    );
    c.font = "10px -apple-system,sans-serif";
  }
  readCanvas.dataset.viewStart = readView.start;
  readCanvas.dataset.viewEnd = readView.end;
  readCanvas.dataset.firstRow = first;
}
$("expand").onclick = () => {
  if (selection) showReads(selection.fragments);
};
$("allreads").onclick = () => showReadIDs(result.reads.map((r) => r.id));
$("close-reads").onclick = () => {
  $("evidence").hidden = true;
  readsOpen = false;
};
startWorker();
draw();
updateButtons();
