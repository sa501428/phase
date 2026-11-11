import { VariantGraphView } from "./variant-graph-view.js";
import { watchVaf } from "../src/live-vaf.js";
import { hiddenBase, hiddenEvent } from "../src/vaf-filter.js";
import {
  IGV_COLORS,
  IGV_ALIGNMENT,
  drawIGVBase,
  drawIGVInsertion,
  drawIGVGap,
} from "./igv-render.js";
import { ungzip } from "pako";
import hg19 from "../references/hg19.js";
import hg38 from "../references/hg38.js";
import { ManualWorkspace } from "./model.js";
import { packPileup, pileupHit } from "./layout.js";
import { groupStatistics, manualHierarchy, siteCounts } from "./statistics.js";
import { hasSequence } from "../src/alignment-geometry.js";
const $ = (id) => document.getElementById(id),
  el = (tag, text, cls) => {
    const e = document.createElement(tag);
    if (text !== undefined) e.textContent = text;
    if (cls) e.className = cls;
    return e;
  };
const colors = IGV_COLORS,
  palette = [
    "#277e91",
    "#8b64ac",
    "#c18a44",
    "#498aa4",
    "#b66077",
    "#698e59",
    "#6875ad",
    "#bc785c",
  ];
const color = (id) => palette[(id - 1) % palette.length],
  pct = (a, b) => (b ? ((100 * a) / b).toFixed(1) + "%" : "unknown");
let worker,
  url,
  timer,
  busy = false,
  opened = false,
  model = null,
  rows = [],
  coverage = null,
  view = null,
  site = null,
  anchor = null,
  inspected = null,
  treeData = null,
  treeHits = [],
  hits = [],
  catalog = null,
  catalogBuild = null,
  genes = [],
  geneModel = null,
  activeGene = null,
  contigs = [],
  lastTooltip = "";
let packed = {
  placements: [],
  byItem: new Map(),
  order: [],
  sections: [],
  height: 570,
  laneCount: 0,
};
let hoverDetails = false,
  contextPoint = null,
  boxSelection = null;
const readHeight = () => +$("readSize").value;
const viewportHeight = () => Math.max(160, $("readScroll").clientHeight || 570);
const pile = $("pileup"),
  track = $("referenceTrack");
const leftMargin = () => (pile.clientWidth < 570 ? 88 : 112);
function status(text, error = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", error);
}
function buttons() {
  const selected = model?.selection.size || 0;
  for (const id of ["merge", "hide", "ungroup", "clearSelection"])
    $(id).disabled = busy || !selected;
  for (const id of ["save", "fasta", "import", "apply", "selectVisible"])
    $(id).disabled = busy || !model;
  $("undo").disabled = busy || !model?.undoStack.length;
  $("redo").disabled = busy || !model?.redoStack.length;
  $("restore").disabled = busy || !model?.hidden.size;
  $("open").disabled = busy;
  $("demo").disabled = busy;
  $("load").disabled = busy || !opened;
  $("selected").textContent =
    `${selected} pair${selected === 1 ? "" : "s"} selected`;
}
function destroyWorker() {
  clearTimeout(timer);
  worker?.terminate();
  if (url) URL.revokeObjectURL(url);
  worker = null;
  url = null;
}
function startWorker() {
  destroyWorker();
  url = URL.createObjectURL(
    new Blob([WORKER_SOURCE], { type: "text/javascript" }),
  );
  worker = new Worker(url);
  worker.onmessage = receive;
  worker.onerror = (e) => {
    clearTimeout(timer);
    busy = false;
    status(e.message, true);
    buttons();
  };
}
function post(message) {
  if (busy) return;
  busy = true;
  buttons();
  timer = setTimeout(() => {
    destroyWorker();
    opened = false;
    busy = false;
    status(
      "Local loading exceeded 90 seconds. Reopen the files and choose a smaller region. No reads were downsampled.",
      true,
    );
    buttons();
  }, 90000);
  worker.postMessage(message);
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
    status(
      `Opened indexed BAM with ${contigs.length} contigs. Enter a gene, exon or coordinates to load reads.`,
    );
  }
  if (m.type === "result") {
    model = new ManualWorkspace(m.data);
    activeGene =
      geneModel &&
      geneModel.chr.replace(/^chr/, "") ===
        model.region.chr.replace(/^chr/, "") &&
      geneModel.exons.some(
        ([s, e]) => s < model.region.end && e > model.region.start,
      )
        ? geneModel
        : null;
    site = null;
    anchor = null;
    inspected = null;
    $("readScroll").scrollTop = 0;
    $("hapName").value = "";
    $("region").value =
      `${model.region.chr}:${model.region.start + 1}-${model.region.end}`;
    $("target").textContent =
      `${model.region.chr}:${(model.region.start + 1).toLocaleString()}–${model.region.end.toLocaleString()}`;
    const first = model.fragments
      .filter(
        (f) =>
          Number.isFinite(f.start) &&
          f.start < model.region.end &&
          f.end > model.region.start,
      )
      .sort((a, b) => a.start - b.start)[0];
    const span = Math.min(600, model.region.end - model.region.start);
    view = {
      start: Math.max(
        model.region.start,
        Math.min(
          (first?.start ?? model.region.start) - 20,
          model.region.end - span,
        ),
      ),
      end: 0,
    };
    view.end = view.start + span;
    $("viewRange").value =
      `${model.region.chr}:${Math.floor(view.start) + 1}-${Math.ceil(view.end)}`;
    model.setSettings(settings());
    document.body.classList.remove("inputs-open");
    document.body.classList.add("top-compact");
    syncChrome();
    refresh();
    inspect(null);
    status(
      `Loaded ${model.fragments.length} read pairs/fragments (${model.reads.length} alignment records). No automatic merging or downsampling.`,
    );
  }
  buttons();
}
function parseRegion(text) {
  const m = text
    .trim()
    .replace(/,/g, "")
    .match(/^([^:\s]+):(\d+)\s*-\s*(\d+)$/);
  if (
    !m ||
    !Number.isSafeInteger(+m[2]) ||
    !Number.isSafeInteger(+m[3]) ||
    +m[2] < 1 ||
    +m[3] < +m[2]
  )
    throw Error("Use contig:start-end with 1-based inclusive coordinates");
  return { chr: m[1], start: +m[2] - 1, end: +m[3] };
}
$("open").onclick = () => {
  try {
    const files = Object.fromEntries(
      ["bam", "index", "reference", "fai"].map((id) => [id, $(id).files[0]]),
    );
    if (!files.bam || !files.index || !files.reference)
      throw Error("Choose BAM, BAI/CSI and matching reference");
    startWorker();
    opened = false;
    clearWorkspace();
    status("Opening local files…");
    post({ type: "open", files });
  } catch (e) {
    status(e.message, true);
  }
};
$("load").onclick = () => {
  try {
    const region = parseRegion($("region").value);
    status("Loading the selected locus and indexed mates…");
    post({ type: "load", region });
  } catch (e) {
    status(e.message, true);
  }
};
$("demo").onclick = () => {
  if (!worker) startWorker();
  geneModel = null;
  status("Preparing synthetic paired reads…");
  post({ type: "demo", scenario: $("scenario").value });
};
function settings() {
  const s = {
    minVaf: +$("minVaf").value / 100,
    mapq: +$("mapq").value,
    baseq: +$("baseq").value,
    excludeHidden: $("excludeHidden").checked,
    includeDuplicates: $("includeDuplicates").checked,
  };
  if (
    !Number.isFinite(s.minVaf) ||
    s.minVaf < 0 ||
    s.minVaf > 1 ||
    !Number.isInteger(s.mapq) ||
    s.mapq < 0 ||
    s.mapq > 254 ||
    !Number.isInteger(s.baseq) ||
    s.baseq < 0 ||
    s.baseq > 93
  )
    throw Error("Check minimum VAF (0–100%), MAPQ and base-quality thresholds");
  return s;
}
function run(action, message) {
  if (!model || busy) return;
  try {
    const value = action();
    anchor = null;
    refresh();
    if (message) status(message);
    return value;
  } catch (e) {
    status(e.message, true);
    buttons();
  }
}
$("apply").onclick = () =>
  run(() => {
    model.setSettings(settings());
    inspect(inspected);
  }, "Evidence settings applied. Existing assignments retained; conflicting consensus sites are flagged.");
watchVaf($("minVaf"), {
  current: () => model?.settings.minVaf,
  busy: () => busy,
  update: (minVaf) =>
    run(() => {
      model.setSettings({ ...model.settings, minVaf });
      $("tooltip").hidden = true;
      inspect(inspected);
    }, "Variant cutoff updated. Existing assignments retained."),
});
$("merge").onclick = () => {
  const id = run(() => {
    const id = model.merge($("hapName").value);
    model.selection.clear();
    model.expanded.delete(id);
    $("hapName").value = "";
    inspected = { kind: "group", id };
    return id;
  }, "Selected pairs merged into a manual haplotype. All pairs deselected. Expand its row to inspect the original read pairs.");
  if (id !== undefined) {
    const index = rows.findIndex(
      (row) => row.kind === "group" && row.group.id === id,
    );
    const p = packed.byItem.get(index);
    if (p) $("readScroll").scrollTop = Math.max(0, p.y - 8);
    drawReads();
    renderAccessible();
    const g = model.groups.find((g) => g.id === id);
    $("mergeFeedback").textContent =
      `Merged ${g.name}: ${g.members.length} pairs. Selection cleared.`;
  } else $("mergeFeedback").textContent = $("status").textContent;
};
$("hide").onclick = () =>
  run(() => {
    model.hideSelection();
    inspected = null;
  }, "Selected pairs hidden. Restore hidden or Undo to bring them back; the evidence checkbox controls whether they contribute to consensus and VAF.");
$("restore").onclick = () =>
  run(() => model.restoreAll(), "Hidden pairs restored.");
$("ungroup").onclick = () =>
  run(() => {
    model.detach();
    inspected = null;
  }, "Selected pairs returned to individual rows.");
$("undo").onclick = () =>
  run(() => {
    model.undo();
    syncSettings();
    inspected = null;
  }, "Last edit undone.");
$("redo").onclick = () =>
  run(() => {
    model.redo();
    syncSettings();
    inspected = null;
  }, "Edit redone.");
function syncSettings() {
  if (!model) return;
  $("minVaf").value = (model.settings.minVaf ?? 0.02) * 100;
  $("mapq").value = model.settings.mapq;
  $("baseq").value = model.settings.baseq;
  $("excludeHidden").checked = model.settings.excludeHidden;
  $("includeDuplicates").checked = model.settings.includeDuplicates;
}
$("clearSelection").onclick = () => {
  if (!model) return;
  model.selection.clear();
  anchor = null;
  drawReads();
  buttons();
  renderAccessible();
};
$("showHidden").onchange = () => {
  if (model) {
    rows = model.rows($("sort").value, site, $("showHidden").checked);
    anchor = null;
    sizeRows();
    drawReads();
    renderAccessible();
  }
};
$("sort").onchange = () => {
  if (!model) return;
  if ($("sort").value === "allele" && site === null) {
    status(
      "Click a reference base or coordinate first to choose the allele-sort site.",
      true,
    );
    $("sort").value = "position";
  }
  rows = model.rows($("sort").value, site, $("showHidden").checked);
  anchor = null;
  $("readScroll").scrollTop = 0;
  sizeRows();
  drawReads();
  renderAccessible();
};
function download(text, name, type = "text/plain") {
  const u = URL.createObjectURL(new Blob([text], { type })),
    a = el("a");
  a.href = u;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(u), 1000);
}
$("save").onclick = () => {
  if (model)
    download(
      model.exportAssignments(),
      "phase-manual-assignments.json",
      "application/json",
    );
};
$("fasta").onclick = () => {
  try {
    download(model.exportFasta(), "phase-manual-haplotypes.fa");
  } catch (e) {
    status(e.message, true);
  }
};
$("import").onchange = async () => {
  const f = $("import").files[0];
  if (!model || !f) return;
  try {
    if (f.size > 16000000) throw Error("Assignment file exceeds 16 MB");
    model.importAssignments(await f.text());
    syncSettings();
    inspected = null;
    refresh();
    status("Manual haplotypes, hidden pairs and evidence settings restored.");
  } catch (e) {
    status(e.message, true);
  } finally {
    $("import").value = "";
  }
};

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
    const name = $("gene").value.trim().toUpperCase();
    if (!name) throw Error("Enter a gene symbol");
    genes = loadCatalog().filter((g) => g[0].toUpperCase() === name);
    if (!genes.length)
      throw Error("Gene not found in the selected offline catalog");
    $("genes").replaceChildren(
      ...genes.map((g, i) => {
        const o = el("option", `${g[0]} · ${g[1]}`);
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
      const o = el("option", t[0]);
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
  geneModel = { name: g[0], chr: g[1], strand: g[2], exons: t[3] };
  const option = el("option", "Whole transcript");
  option.value = "";
  $("exon").replaceChildren(
    option,
    ...t[3].map(([s, e], i) => {
      const o = el(
        "option",
        `Exon ${g[2] === "-" ? t[3].length - i : i + 1} · ${s + 1}–${e}`,
      );
      o.value = i;
      return o;
    }),
  );
  setGeneRegion([t[1], t[2]]);
}
function setGeneRegion([s, e]) {
  const g = genes[+$("genes").value],
    aliases = [g[1], g[1].startsWith("chr") ? g[1].slice(3) : "chr" + g[1]],
    chr = contigs.find((c) => aliases.includes(c.refName))?.refName || g[1];
  $("region").value = `${chr}:${s + 1}-${e}`;
  status(
    `Selected ${g[0]} (${$("build").value}). Load the selected region to begin a new workspace.`,
  );
}
$("genes").onchange = chooseGene;
$("transcript").onchange = chooseTranscript;
$("exon").onchange = () => {
  const t = genes[+$("genes").value]?.[6][+$("transcript").value];
  if (t)
    setGeneRegion(
      $("exon").value === "" ? [t[1], t[2]] : t[3][+$("exon").value],
    );
};
$("build").onchange = () => {
  catalog = null;
  catalogBuild = null;
  genes = [];
  geneModel = null;
  for (const [id, text] of [
    ["genes", "Gene locus"],
    ["transcript", "Transcript"],
    ["exon", "Whole transcript"],
  ]) {
    const o = el("option", text);
    o.value = "";
    $(id).replaceChildren(o);
  }
};
function clearWorkspace() {
  model = null;
  variantNetwork.update(null);
  activeGene = null;
  packed = {
    placements: [],
    byItem: new Map(),
    order: [],
    sections: [],
    height: 570,
    laneCount: 0,
  };
  hoverDetails = false;
  contextPoint = null;
  boxSelection = null;
  $("readMenu").hidden = true;
  $("readMenu").replaceChildren();
  rows = [];
  coverage = null;
  view = null;
  site = null;
  anchor = null;
  inspected = null;
  treeData = null;
  treeHits = [];
  hits = [];
  lastTooltip = "";
  $("tooltip").hidden = true;
  $("tooltip").textContent = "";
  $("tree").title = "";
  $("siteStats").replaceChildren();
  $("siteTitle").textContent = "";
  $("hapName").value = "";
  $("mergeFeedback").textContent = "";
  for (const key of ["firstRow", "viewStart", "viewEnd", "rowCount"])
    delete pile.dataset[key];
  $("readScroll").scrollTop = 0;
  sizeRows();
  $("groupList").replaceChildren();
  $("inspector").replaceChildren(
    el("p", "Select a read pair or constructed haplotype.", "muted"),
  );
  $("inspectorTitle").textContent = "Select a pair or haplotype";
  $("target").textContent = "Manual haplotype workspace";
  $("counts").textContent = "Pairs stay separate until you merge them.";
  $("treeSummary").textContent =
    "Merge read pairs to build groups for comparison.";
  $("sitePanel").hidden = true;
  $("siteLabel").textContent = "No site selected";
  $("viewRange").value = "";
  $("rowCount").textContent = "0 rows";
  $("accessibleRows").replaceChildren();
  drawReads();
  drawTrack();
  drawTree();
  buttons();
}
$("reset").onclick = () => {
  destroyWorker();
  busy = false;
  opened = false;
  contigs = [];
  catalog = null;
  catalogBuild = null;
  genes = [];
  geneModel = null;
  for (const id of [
    "bam",
    "index",
    "reference",
    "fai",
    "import",
    "region",
    "gene",
    "hapName",
  ])
    $(id).value = "";
  $("showHidden").checked = false;
  $("includeDuplicates").checked = false;
  document.body.classList.add("inputs-open");
  document.body.classList.remove("top-compact");
  $("sort").value = "position";
  for (const [id, text] of [
    ["genes", "Gene locus"],
    ["transcript", "Transcript"],
    ["exon", "Whole transcript"],
  ]) {
    const o = el("option", text);
    o.value = "";
    $(id).replaceChildren(o);
  }
  clearWorkspace();
  syncChrome();
  startWorker();
  status("Session cleared. No case data was saved in browser storage.");
};
function syncChrome() {
  $("toggleInputs").setAttribute(
    "aria-expanded",
    String(document.body.classList.contains("inputs-open")),
  );
  $("toggleAnalysis").setAttribute(
    "aria-expanded",
    String(document.body.classList.contains("analysis-open")),
  );
  const expanded = !document.body.classList.contains("top-compact");
  $("toggleTop").setAttribute("aria-expanded", String(expanded));
  $("toggleTop").textContent = expanded ? "Compact top" : "Expand top";
  requestAnimationFrame(() => {
    sizeRows();
    drawTrack();
    drawReads();
    drawTree();
    renderAccessible();
  });
}
$("toggleInputs").onclick = () => {
  document.body.classList.toggle("inputs-open");
  syncChrome();
};
$("toggleAnalysis").onclick = () => {
  document.body.classList.toggle("analysis-open");
  syncChrome();
};
$("toggleTop").onclick = () => {
  document.body.classList.toggle("top-compact");
  syncChrome();
};
$("readSize").onchange = () => {
  sizeRows();
  drawReads();
  renderAccessible();
};
$("includeDuplicates").onchange = () => {
  if (!model) return;
  model.setSettings(settings());
  if (!model.settings.includeDuplicates)
    for (const id of model.selection)
      if (model.fragments[id].duplicate) model.selection.delete(id);
  anchor = null;
  refresh();
};
function refresh() {
  if (model) pile.dataset.minVaf = model.settings.minVaf;
  if (!model) return;
  rows = model.rows($("sort").value, site, $("showHidden").checked);
  coverage = model.coverage();
  treeData = manualHierarchy(model, +$("treeCut").value / 100);
  sizeRows();
  $("counts").textContent =
    `${model.fragments.length} pairs/fragments · ${model.reads.length} records · ${model.groups.length} manual haplotypes · ${model.hidden.size} hidden pairs`;

  drawReads();
  drawTrack();
  renderGroups();
  drawTree();
  inspect(inspected);
  renderSite();
  renderAccessible();
  buttons();
  variantNetwork.update(model);
}
function sizeRows() {
  const height = viewportHeight();
  if (model && view)
    packed = packPileup(model, rows, {
      view,
      width: pile.clientWidth || 700,
      left: leftMargin(),
      readHeight: readHeight(),
      sort: $("sort").value,
      site,
      minHeight: height,
    });
  const rule = [...document.styleSheets[0].cssRules].find(
    (r) => r.selectorText === "#readSpace",
  );
  if (rule) rule.style.height = Math.max(height, packed.height) + "px";
  pile.style.height = height + "px";
  if (model) {
    const visible = packed.placements.filter((p) => !p.header).length,
      dups = model.fragments.filter((f) => f.duplicate).length;
    $("rowCount").textContent =
      `${packed.laneCount} packed lanes · ${visible} visible pairs${!model.settings.includeDuplicates && dups ? ` · ${dups} duplicate-flagged pairs filtered` : ""}`;
    $("compactTarget").textContent =
      `${model.region.chr}:${model.region.start + 1}–${model.region.end} · ${model.groups.length} haplotypes`;
  }
}
function resize(canvas, height) {
  const w = canvas.clientWidth || 700,
    d = devicePixelRatio || 1;
  canvas.width = Math.round(w * d);
  canvas.height = Math.round(height * d);
  const c = canvas.getContext("2d");
  c.setTransform(d, 0, 0, d, 0, 0);
  return { c, w, h: height };
}
const variantNetwork = new VariantGraphView($("variantGraph"), {
  selectPairs(ids) {
    if (!model) return;
    model.selection = new Set(model.activeIDs(ids));
    anchor = null;
    drawReads();
    buttons();
    renderAccessible();
    status(
      `Selected ${model.selection.size} supporting pairs from the variant graph.`,
    );
  },
  focusVariant(v) {
    if (!model || !v) return;
    site = v.pos;
    const span = Math.min(
      model.region.end - model.region.start,
      Math.max(120, view.end - view.start),
    );
    setView(v.pos - span / 2, span);
    renderSite();
  },
});
const axisWidth = () =>
  Math.min(track.clientWidth || 700, pile.clientWidth || 700);
const xPosition = (p, w) =>
  leftMargin() +
  ((p - view.start) / (view.end - view.start)) *
    (axisWidth() - leftMargin() - 12);
function setView(start, span) {
  if (!model) return;
  const full = model.region.end - model.region.start;
  span = Math.max(Math.min(12, full), Math.min(full, span));
  start = Math.max(
    model.region.start,
    Math.min(model.region.end - span, start),
  );
  view = { start, end: start + span };
  closeReadMenu();
  sizeRows();
  $("viewRange").value =
    `${model.region.chr}:${Math.floor(start) + 1}-${Math.ceil(view.end)}`;
  drawTrack();
  drawReads();
}
function zoom(f, center = 0.5) {
  if (!model) return;
  const span = view.end - view.start,
    next = span * f;
  setView(view.start + span * center - next * center, next);
}
$("left").onclick = () => {
  if (model)
    setView(view.start - (view.end - view.start) * 0.3, view.end - view.start);
};
$("right").onclick = () => {
  if (model)
    setView(view.start + (view.end - view.start) * 0.3, view.end - view.start);
};
$("zoomIn").onclick = () => zoom(0.5);
$("zoomOut").onclick = () => zoom(2);
$("fit").onclick = () => {
  if (model) setView(model.region.start, model.region.end - model.region.start);
};
$("go").onclick = () => {
  try {
    const r = parseRegion($("viewRange").value);
    if (
      !model ||
      r.chr !== model.region.chr ||
      r.start < model.region.start ||
      r.end > model.region.end
    )
      throw Error(
        "View must stay inside the loaded target. Enter new coordinates at left and load to change the workspace.",
      );
    setView(r.start, r.end - r.start);
  } catch (e) {
    status(e.message, true);
  }
};
function drawTrack() {
  const { c, w } = resize(track, 128),
    left = leftMargin();
  c.fillStyle = "#fcfdfe";
  c.fillRect(0, 0, w, 128);
  c.fillStyle = "#8196a4";
  c.font = "10px -apple-system,sans-serif";
  c.fillText("Coordinates", 12, 17);
  c.fillText("Pair coverage", 12, 50);
  c.fillText("Reference", 12, 89);
  c.fillText(activeGene?.name || "Transcript", 12, 118);
  if (!model || !view) return;
  const span = view.end - view.start,
    ppb = (axisWidth() - left - 12) / span,
    step = Math.max(1, Math.pow(10, Math.floor(Math.log10(span / 5)))),
    tick = span / step > 12 ? step * 5 : step;
  c.save();
  c.beginPath();
  c.rect(left, 0, w - left, 128);
  c.clip();
  for (let p = Math.ceil(view.start / tick) * tick; p < view.end; p += tick) {
    const x = xPosition(p, w);
    c.fillStyle = "#6d889b";
    c.fillText((p + 1).toLocaleString(), x + 3, 17);
    c.strokeStyle = "#e5edf3";
    c.beginPath();
    c.moveTo(x, 21);
    c.lineTo(x, 104);
    c.stroke();
  }
  let max = 0;
  const start = Math.max(0, Math.floor(view.start - model.region.start)),
    end = Math.min(coverage.length, Math.ceil(view.end - model.region.start));
  for (let i = start; i < end; i++) max = Math.max(max, coverage[i]);
  for (let x = left; x < w - 12; x++) {
    const a = Math.max(
        0,
        Math.floor(view.start - model.region.start + (x - left) / ppb),
      ),
      b = Math.min(
        coverage.length,
        Math.ceil(view.start - model.region.start + (x - left + 1) / ppb),
      );
    let n = 0;
    for (let p = a; p < b; p++) n = Math.max(n, coverage[p]);
    if (n) {
      c.fillStyle = "#9fbccc";
      c.fillRect(
        x,
        67 - (35 * n) / Math.max(1, max),
        1,
        (35 * n) / Math.max(1, max),
      );
    }
  }
  for (let p = Math.floor(view.start); p < view.end; p++) {
    const b = model.reference[p - model.region.start],
      x = xPosition(p, w);
    if (!b) continue;
    drawIGVBase(c, b, x, 76, ppb, 18);
  }
  if (
    activeGene &&
    activeGene.chr.replace(/^chr/, "") === model.region.chr.replace(/^chr/, "")
  ) {
    c.strokeStyle = "#95aebd";
    c.beginPath();
    c.moveTo(left, 115);
    c.lineTo(w - 12, 115);
    c.stroke();
    for (const [s, e] of activeGene.exons) {
      c.fillStyle = "#4b839c";
      c.fillRect(
        xPosition(s, w),
        109,
        Math.max(2, xPosition(e, w) - xPosition(s, w)),
        12,
      );
    }
  }
  if (site !== null) {
    c.fillStyle = "#e5b64744";
    c.fillRect(xPosition(site, w), 0, Math.max(2, ppb), 105);
  }
  c.restore();
  c.fillStyle = "#8097a5";
  c.font = "9px -apple-system,sans-serif";
  c.fillText(`${max}× max`, 12, 65);
}
function drawReads() {
  const { c, w } = resize(pile, viewportHeight()),
    height = viewportHeight(),
    left = leftMargin(),
    scroll = $("readScroll").scrollTop;
  c.fillStyle = "#ffffff";
  c.fillRect(0, 0, w, height);
  hits = [];
  if (!model || !view) {
    c.fillStyle = "#8797a1";
    c.font = "13px -apple-system,sans-serif";
    c.textAlign = "center";
    c.fillText(
      "Open files and choose a gene or interval to load its pileup",
      w / 2,
      Math.min(160, height / 2),
    );
    return;
  }
  pile.dataset.firstRow = Math.floor(scroll / (readHeight() + 9));
  pile.dataset.viewStart = view.start;
  pile.dataset.viewEnd = view.end;
  pile.dataset.rowCount = rows.length;
  pile.dataset.laneCount = packed.laneCount;
  pile.dataset.readHeight = readHeight();
  pile.dataset.leftMargin = left;
  const ppb = (w - left - 12) / (view.end - view.start);
  c.fillStyle = "#f8fafb";
  c.fillRect(0, 0, left - 1, height);
  if (site !== null) {
    c.fillStyle = "#f0cf7733";
    c.fillRect(xPosition(site, w), 0, Math.max(2, ppb), height);
  }
  for (const section of packed.sections)
    if (section.end >= scroll && section.y < scroll + height) {
      c.fillStyle = "#8295a1";
      c.font = "10px -apple-system,sans-serif";
      c.textAlign = "left";
      c.fillText(section.label, 8, Math.max(14, section.y - scroll + 14));
    }
  for (const p of packed.placements) {
    if (p.y + p.height < scroll || p.y > scroll + height) continue;
    const row = rows[p.index],
      y = p.y - scroll,
      selected = row.members.some((id) => model.selection.has(id));
    if (p.header) {
      c.fillStyle = "#f0f5f8";
      c.fillRect(0, y, w, p.height);
      c.fillStyle = color(row.group.id);
      c.font = "bold 11px -apple-system,sans-serif";
      c.textAlign = "left";
      c.fillText(
        (model.expanded.has(row.group.id) ? "▾ " : "▸ ") +
          row.group.name.slice(0, 12),
        7,
        y + 13,
      );
      c.fillStyle = "#8397a4";
      c.font = "9px -apple-system,sans-serif";
      c.fillText(
        `${model.activeIDs(row.group.members).length} pairs`,
        8,
        y + 26,
      );
    }
    c.save();
    c.beginPath();
    c.rect(left, y, w - left, p.height + 2);
    c.clip();
    if (row.kind !== "group" && model.hidden.has(row.fragment.id))
      c.globalAlpha = 0.32;
    if (selected) {
      c.fillStyle = "#d8eaff";
      c.fillRect(p.x1 - 3, y - 1, Math.max(10, p.x2 - p.x1 + 6), p.height + 2);
    }
    if (row.kind === "group") drawConsensus(c, row.group, y, w, ppb);
    else drawPair(c, row.fragment, y, w, ppb);
    if (selected) {
      c.strokeStyle = "#1d68a5";
      c.lineWidth = 1.4;
      c.strokeRect(p.x1 - 3, y, Math.max(10, p.x2 - p.x1 + 6), p.height);
    }
    c.restore();
  }
  if (!packed.placements.length) {
    c.fillStyle = "#8697a2";
    c.textAlign = "center";
    c.font = "12px sans-serif";
    c.fillText(
      model.hidden.size
        ? "No displayed pairs here. Restore hidden pairs or pan to another position."
        : "No displayed pairs in this view. Pan or change duplicate filtering.",
      (w + left) / 2,
      100,
    );
  }
  c.strokeStyle = "#dbe3e8";
  c.beginPath();
  c.moveTo(left - 0.5, 0);
  c.lineTo(left - 0.5, height);
  c.stroke();
  if (boxSelection) {
    c.fillStyle = "#3e88c228";
    c.strokeStyle = "#3e88c2";
    const { x, y, endX, endY } = boxSelection;
    c.fillRect(
      Math.min(x, endX),
      Math.min(y, endY),
      Math.abs(endX - x),
      Math.abs(endY - y),
    );
    c.strokeRect(
      Math.min(x, endX),
      Math.min(y, endY),
      Math.abs(endX - x),
      Math.abs(endY - y),
    );
  }
}
function drawPair(c, f, y, w, ppb) {
  const rs = f.readIds.map((id) => model.readById.get(id)),
    mapped = rs.filter((r) => r.alignmentUsable),
    h = readHeight(),
    yy = y,
    center = yy + h / 2;
  if (mapped.length > 1) {
    const x1 = xPosition(Math.min(...mapped.map((r) => r.start)), w),
      x2 = xPosition(Math.max(...mapped.map((r) => r.end)), w);
    c.strokeStyle = IGV_ALIGNMENT;
    c.lineWidth = 1;
    c.beginPath();
    c.moveTo(x1, center);
    c.lineTo(x2, center);
    c.stroke();
    hits.push({ fragment: f.id, x: x1, y: yy, w: x2 - x1, h });
  }
  for (const r of rs) {
    c.save();
    if (r.flags & (256 | 512)) c.globalAlpha *= 0.45;
    if (!r.alignmentUsable) {
      const x = xPosition(Math.max(model.region.start, r.start), w);
      c.fillStyle = "#b3b8be";
      c.fillRect(x, yy, 9, h);
      c.fillStyle = "#fff";
      c.font = "11px sans-serif";
      c.fillText("?", x + 1, yy + h - 3);
      hits.push({ read: r.id, fragment: f.id, x, y: yy, w: 10, h });
      c.restore();
      continue;
    }
    for (const run of r.runs) {
      const { op, length: n, refStart: p, queryStart: q } = run,
        x = xPosition(p, w),
        end = xPosition(p + n, w);
      if ("M=X".includes(op)) {
        c.fillStyle = IGV_ALIGNMENT;
        if (hasSequence(r) && r.mq !== 0)
          c.fillRect(x, yy, Math.max(1, end - x), h);
        else {
          c.strokeStyle = "#a8adb2";
          c.strokeRect(x, yy, Math.max(1, end - x), h);
        }
        hits.push({
          read: r.id,
          fragment: f.id,
          x,
          y: yy,
          w: Math.max(3, end - x),
          h,
        });
        if (hasSequence(r))
          for (let j = 0; j < n; j++) {
            const pos = p + j,
              b = r.seq[q + j],
              ref = model.reference[pos - model.region.start],
              xx = xPosition(pos, w);
            if (xx < leftMargin() - ppb || xx > w) continue;
            const mismatch = ref && b !== ref && b !== "N";
            if (mismatch && !hiddenBase(r, q + j))
              drawIGVBase(c, b, xx, yy, ppb, h, r.qual?.[q + j]);
          }
      } else if (op === "I" && !hiddenEvent(r, op, p)) {
        const glyph = drawIGVInsertion(c, x, yy, h, n, ppb);
        hits.push({
          read: r.id,
          fragment: f.id,
          kind: "insertion",
          x: glyph.left - 2,
          y: yy,
          w: glyph.right - glyph.left + 4,
          h,
          detail: `Insertion at ${p + 1}: ${hasSequence(r) ? r.seq.slice(q, q + n) : "sequence unavailable"}`,
        });
      } else if ((op === "D" || op === "N") && !hiddenEvent(r, op, p)) {
        drawIGVGap(c, x, end, yy, h, n, op);
        hits.push({
          read: r.id,
          fragment: f.id,
          x,
          y: yy,
          w: Math.max(3, end - x),
          h,
          detail: `${n}${op} at ${p + 1}`,
        });
      } else if (op === "S" && hasSequence(r)) {
        const side = q === 0 ? "left" : "right",
          clip = { pos: p, side, seq: r.seq.slice(q, q + n) };
        drawClip(c, clip, yy, w, ppb, r.qual?.slice(q, q + n));
        const from = side === "left" ? p - n : p;
        hits.push({
          read: r.id,
          fragment: f.id,
          x: xPosition(from, w),
          y: yy,
          w: Math.max(6, n * ppb),
          h,
          detail: `${side === "left" ? "Leading" : "Trailing"} soft clip (${n} bp, unaligned): ${clip.seq}`,
        });
      }
    }
    c.fillStyle = IGV_ALIGNMENT;
    c.lineWidth = 1;
    const arrow = xPosition(r.flags & 16 ? r.start : r.end, w),
      dir = r.flags & 16 ? -1 : 1;
    c.beginPath();
    c.moveTo(arrow, yy);
    c.lineTo(arrow + dir * 4, center);
    c.lineTo(arrow, yy + h);
    c.fill();
    c.restore();
  }
}
function drawClip(c, clip, yy, w, ppb, qualities) {
  const h = readHeight(),
    from = clip.side === "left" ? clip.pos - clip.seq.length : clip.pos;
  c.fillStyle = IGV_ALIGNMENT;
  c.fillRect(xPosition(from, w), yy, clip.seq.length * ppb, h);
  for (let j = 0; j < clip.seq.length; j++) {
    const x = xPosition(from + j, w);
    if (x + ppb < leftMargin() || x > w) continue;
    drawIGVBase(c, clip.seq[j], x, yy, ppb, h, qualities?.[j]);
  }
}
function drawConsensus(c, g, y, w, ppb) {
  const s = model.groupSummary(g),
    h = readHeight(),
    yy = y + 5;
  // Fill covered runs once so fractional per-base edges cannot produce a
  // barcode appearance at overview zoom. Unknown gaps remain unpainted.
  let runStart = null,
    runEnd = null;
  const fillRun = () => {
    if (runStart !== null) {
      c.fillStyle = color(g.id);
      c.fillRect(
        xPosition(runStart, w),
        yy,
        Math.max(1, xPosition(runEnd, w) - xPosition(runStart, w)),
        h,
      );
      runStart = null;
    }
  };
  for (const call of s.bases) {
    const deletion = !call.conflict && call.values[0].value === "-";
    if (deletion) {
      fillRun();
      continue;
    }
    if (runStart === null) {
      runStart = call.pos;
      runEnd = call.pos + 1;
    } else if (call.pos === runEnd) runEnd++;
    else {
      fillRun();
      runStart = call.pos;
      runEnd = call.pos + 1;
    }
  }
  fillRun();
  for (const call of s.bases) {
    const p = call.pos,
      x = xPosition(p, w),
      b = call.conflict ? "N" : call.values[0].value,
      ref = model.reference[p - model.region.start];
    if (x < leftMargin() - ppb || x > w) continue;
    if (b !== "-" && (call.conflict || b !== ref))
      drawIGVBase(c, b, x, yy, ppb, h);
  }
  const deleted = s.bases.filter(
    (call) => !call.conflict && call.values[0].value === "-",
  );
  for (let i = 0; i < deleted.length;) {
    const start = deleted[i].pos;
    let end = start + 1;
    i++;
    while (i < deleted.length && deleted[i].pos === end) {
      end++;
      i++;
    }
    drawIGVGap(
      c,
      xPosition(start, w),
      xPosition(end, w),
      yy,
      h,
      end - start,
      "D",
    );
  }
  for (const ins of s.insertions) {
    const length = Math.max(1, ...ins.values.map((v) => v.value.length)),
      glyph = drawIGVInsertion(c, xPosition(ins.pos, w), yy, h, length, ppb);
    hits.push({
      kind: "insertion",
      group: g.id,
      x: glyph.left - 2,
      y: yy,
      w: glyph.right - glyph.left + 4,
      h,
      detail: `Insertion at ${model.region.chr}:${ins.pos + 1}: ${ins.values.map((v) => v.value || "no insertion").join(" / ")}`,
    });
  }
  for (const clip of s.clips) drawClip(c, clip, yy, w, ppb);
  if (!s.bases.length) {
    c.fillStyle = "#8b9ca7";
    c.font = "11px sans-serif";
    c.fillText(
      "No callable consensus in current evidence",
      leftMargin() + 8,
      yy + h - 3,
    );
  }
  if (s.conflicts.length) {
    c.fillStyle = "#af791f";
    c.font = "10px sans-serif";
    c.fillText(`${s.conflicts.length} conflicts`, w - 90, y + 29);
  }
}
function rowAt(e) {
  const b = pile.getBoundingClientRect(),
    x = e.clientX - b.left,
    y = e.clientY - b.top,
    worldY = y + $("readScroll").scrollTop,
    hit =
      x < leftMargin()
        ? packed.placements.find(
            (p) => p.header && worldY >= p.y && worldY <= p.y + p.height,
          )
        : pileupHit(packed, x, worldY);
  return { index: hit?.index ?? -1, x, y, hit };
}

function chooseRow(index, shift = false) {
  const row = rows[index];
  if (!row) return;
  if (shift && anchor !== null) {
    const order = packed.order;
    const target = order.indexOf(index),
      a = order.findIndex((i) => rows[i].key === anchor);
    for (
      let i = Math.min(a < 0 ? target : a, target);
      i <= Math.max(a < 0 ? target : a, target);
      i++
    )
      for (const id of rows[order[i]].members) model.selection.add(id);
  } else {
    const all =
      row.members.length && row.members.every((id) => model.selection.has(id));
    for (const id of row.members)
      if (all) model.selection.delete(id);
      else model.selection.add(id);
    anchor = row.key;
  }
  inspect(
    row.kind === "group"
      ? { kind: "group", id: row.group.id }
      : { kind: "fragment", id: row.fragment.id },
  );
  drawReads();
  buttons();
  renderAccessible();
}
let drag = null;
pile.onpointerdown = (e) => {
  if (!model || e.button !== 0) return;
  closeReadMenu();
  pile.focus();
  pile.setPointerCapture(e.pointerId);
  const b = pile.getBoundingClientRect();
  drag = {
    x: e.clientX,
    y: e.clientY,
    start: view.start,
    span: view.end - view.start,
    moved: false,
    select: e.shiftKey,
  };
  if (e.shiftKey)
    boxSelection = {
      x: e.clientX - b.left,
      y: e.clientY - b.top,
      endX: e.clientX - b.left,
      endY: e.clientY - b.top,
    };
};
pile.onpointermove = (e) => {
  if (drag) {
    const dx = e.clientX - drag.x,
      dy = e.clientY - drag.y;
    if (Math.abs(dx) + (drag.select ? Math.abs(dy) : 0) > 4) drag.moved = true;
    if (drag.moved) {
      if (drag.select) {
        const b = pile.getBoundingClientRect();
        boxSelection.endX = e.clientX - b.left;
        boxSelection.endY = e.clientY - b.top;
        drawReads();
      } else
        setView(
          drag.start -
            (dx / (pile.clientWidth - leftMargin() - 12)) * drag.span,
          drag.span,
        );
    }
    return;
  }
  if ($("readMenu").hidden) {
    const insertion = insertionAt(e);
    if (insertion) showInsertion(insertion, e);
    else if (hoverDetails) tooltip(e);
    else if (!contextPoint) {
      $("tooltip").hidden = true;
      lastTooltip = "";
    }
  }
};
pile.onpointerup = (e) => {
  if (!drag) return;
  const clickedInsertion = !drag.moved ? insertionAt(e) : null;
  if (drag.select && drag.moved && boxSelection) {
    const b = boxSelection,
      scroll = $("readScroll").scrollTop,
      x1 = Math.min(b.x, b.endX),
      x2 = Math.max(b.x, b.endX),
      y1 = Math.min(b.y, b.endY) + scroll,
      y2 = Math.max(b.y, b.endY) + scroll;
    for (const p of packed.placements)
      if (p.x1 <= x2 && p.x2 >= x1 && p.y <= y2 && p.y + p.height >= y1)
        for (const id of rows[p.index].members) model.selection.add(id);
    anchor = null;
  } else if (!drag.moved) {
    const { index, x } = rowAt(e),
      row = rows[index];
    if (row?.kind === "group" && x < leftMargin()) {
      if (model.expanded.has(row.group.id)) model.expanded.delete(row.group.id);
      else model.expanded.add(row.group.id);
      rows = model.rows($("sort").value, site, $("showHidden").checked);
      anchor = null;
      sizeRows();
      renderAccessible();
    } else if (!clickedInsertion || e.shiftKey) chooseRow(index, e.shiftKey);
  }
  boxSelection = null;
  drag = null;
  drawReads();
  buttons();
  renderAccessible();
  if (clickedInsertion) {
    contextPoint = { clientX: e.clientX, clientY: e.clientY };
    showInsertion(clickedInsertion, e);
  }
};
pile.onpointercancel = () => {
  drag = null;
  boxSelection = null;
  drawReads();
};
pile.onpointerleave = () => {
  if (!contextPoint) $("tooltip").hidden = true;
  lastTooltip = "";
};
pile.onwheel = (e) => {
  if (!model || (!e.ctrlKey && !e.metaKey)) return;
  e.preventDefault();
  const b = pile.getBoundingClientRect();
  zoom(
    Math.exp(e.deltaY * 0.003),
    Math.max(
      0,
      Math.min(
        1,
        (e.clientX - b.left - leftMargin()) / (pile.clientWidth - leftMargin()),
      ),
    ),
  );
};
$("readScroll").onscroll = () => {
  closeReadMenu();
  drawReads();
  renderAccessible();
};
$("selectVisible").onclick = () => {
  if (!model) return;
  const scroll = $("readScroll").scrollTop,
    height = viewportHeight();
  for (const p of packed.placements)
    if (p.y + p.height >= scroll && p.y < scroll + height)
      for (const id of rows[p.index].members) model.selection.add(id);
  drawReads();
  buttons();
  renderAccessible();
};
function closeReadMenu() {
  contextPoint = null;
  $("readMenu").hidden = true;
  $("tooltip").hidden = true;
  lastTooltip = "";
}
function readMenu(e) {
  if (!model) return;
  e.preventDefault();
  drag = null;
  boxSelection = null;
  const { index } = rowAt(e);
  if (index < 0) {
    closeReadMenu();
    return;
  }
  const row = rows[index],
    menu = $("readMenu");
  contextPoint = { clientX: e.clientX, clientY: e.clientY };
  menu.replaceChildren();
  const title = el(
    "div",
    row.kind === "group" ? row.group.name : row.fragment.name,
    "menu-title",
  );
  menu.append(title);
  const option = (name, fn) => {
    const b = el("button", name);
    b.setAttribute("role", "menuitem");
    b.onclick = () => {
      menu.hidden = true;
      fn();
    };
    menu.append(b);
  };
  option(row.kind === "group" ? "Haplotype details" : "Read details", () => {
    tooltip(contextPoint, true);
    inspect(
      row.kind === "group"
        ? { kind: "group", id: row.group.id }
        : { kind: "fragment", id: row.fragment.id },
    );
  });
  option(
    hoverDetails ? "Disable hover details" : "Enable hover details",
    () => {
      hoverDetails = !hoverDetails;
      closeReadMenu();
    },
  );
  option(
    row.members.every((id) => model.selection.has(id))
      ? "Deselect this pair"
      : "Select this pair",
    () => {
      chooseRow(index, false);
      closeReadMenu();
    },
  );
  menu.hidden = false;
  $("tooltip").hidden = true;
  menu.style.left = Math.max(8, Math.min(innerWidth - 245, e.clientX)) + "px";
  menu.style.top = Math.max(8, Math.min(innerHeight - 170, e.clientY)) + "px";
  menu.querySelector("button")?.focus();
}
pile.oncontextmenu = readMenu;
document.addEventListener("pointerdown", (e) => {
  if (!$("readMenu").contains(e.target) && e.target !== pile) closeReadMenu();
});

function insertionAt(e) {
  const b = pile.getBoundingClientRect(),
    x = e.clientX - b.left,
    y = e.clientY - b.top;
  return [...hits]
    .reverse()
    .find(
      (h) =>
        h.kind === "insertion" &&
        x >= h.x &&
        x <= h.x + h.w &&
        y >= h.y &&
        y <= h.y + h.h,
    );
}
function showInsertion(hit, e) {
  $("tooltip").textContent = hit.detail;
  $("tooltip").hidden = false;
  lastTooltip = "insertion:" + hit.detail;
  $("tooltip").style.left =
    Math.max(8, Math.min(innerWidth - 390, e.clientX + 14)) + "px";
  $("tooltip").style.top =
    Math.max(8, Math.min(innerHeight - 180, e.clientY + 15)) + "px";
}
function tooltip(e, force = false) {
  if (!model || (!hoverDetails && !force)) return;
  const { index, x, y } = rowAt(e),
    row = rows[index];
  if (!row) {
    $("tooltip").hidden = true;
    return;
  }
  const pos = Math.floor(
      view.start +
        ((x - leftMargin()) / (pile.clientWidth - leftMargin() - 12)) *
          (view.end - view.start),
    ),
    hit = [...hits]
      .reverse()
      .find((h) => x >= h.x && x <= h.x + h.w && y >= h.y && y <= h.y + h.h),
    key = row.key + ":" + pos + ":" + (hit?.read ?? "");
  let text;
  if (row.kind === "group") {
    const g = row.group,
      s = model.groupSummary(g),
      rs = g.members
        .flatMap((id) => model.fragments[id].readIds)
        .map((id) => model.readById.get(id)),
      mq = rs.filter((r) => r.mq !== 255).map((r) => r.mq);
    text = `${g.name}\n${g.members.length} assigned pairs · ${g.members.filter((id) => model.hidden.has(id)).length} hidden\nMAPQ range: ${mq.length ? mq.reduce((a, b) => Math.min(a, b), Infinity) + "–" + mq.reduce((a, b) => Math.max(a, b), 0) : "unknown"}\n${s.conflicts.length} conflicting sites · ${s.unknown} pairs without callable sequence`;
    if (pos >= model.region.start && pos < model.region.end) {
      const stats = siteCounts(model, g.members, pos);
      text += `\n${model.region.chr}:${pos + 1} · ${[...stats.counts].map(([b, n]) => `${b}: ${n}/${stats.total}`).join(" · ") || "uncovered"}`;
    }
  } else {
    const f = row.fragment,
      r =
        hit?.read !== undefined
          ? model.readById.get(hit.read)
          : model.readById.get(f.readIds[0]);
    text = `${f.name}\n${r.flags & 128 ? "Read 2" : r.flags & 64 ? "Read 1" : "Read"} · MAPQ: ${r.mq === 255 ? "unknown (255)" : r.mq}\n${r.chr || model.region.chr}:${r.start + 1}–${r.end} · ${r.flags & 16 ? "reverse" : "forward"}\nCIGAR: ${r.cigar || "unavailable"} · flags: ${r.flags}${r.flags & 1024 ? " · duplicate" : ""}${r.flags & 2048 ? " · supplementary" : ""}${r.flags & 256 ? " · secondary" : ""}\n${f.readIds
      .map((id) => {
        const r = model.readById.get(id);
        return `${r.flags & 128 ? "R2" : r.flags & 64 ? "R1" : "read"} MAPQ ${r.mq === 255 ? "unknown" : r.mq}`;
      })
      .join(" · ")}`;
    if (r.displayIssue) text += "\n" + r.displayIssue;
    if (r.mateChr) text += `\nMate: ${r.mateChr}:${r.mate + 1}`;
    if (hit?.detail) text += "\n" + hit.detail;
    else if (pos >= model.region.start && pos < model.region.end) {
      for (const run of r.runs)
        if (
          "M=X".includes(run.op) &&
          pos >= run.refStart &&
          pos < run.refStart + run.length
        ) {
          const q = run.queryStart + pos - run.refStart;
          text += `\n${pos + 1}: ${hasSequence(r) ? r.seq[q] : "SEQ unavailable"} · reference ${model.reference[pos - model.region.start]} · base Q ${r.qual?.[q] === 255 ? "unknown" : (r.qual?.[q] ?? "unknown")}`;
        }
    }
  }
  if (key !== lastTooltip) {
    $("tooltip").textContent = text;
    lastTooltip = key;
  }
  $("tooltip").hidden = false;
  $("tooltip").style.left =
    Math.min(innerWidth - 390, Math.max(8, e.clientX + 14)) + "px";
  $("tooltip").style.top =
    Math.min(innerHeight - 240, Math.max(8, e.clientY + 15)) + "px";
}
track.onclick = (e) => {
  if (!model) return;
  const b = track.getBoundingClientRect(),
    x = e.clientX - b.left;
  if (x < leftMargin() || x >= axisWidth() - 12) return;
  site = Math.floor(
    view.start +
      ((x - leftMargin()) / (axisWidth() - leftMargin() - 12)) *
        (view.end - view.start),
  );
  site = Math.max(model.region.start, Math.min(model.region.end - 1, site));
  $("siteLabel").textContent = `${model.region.chr}:${site + 1}`;
  drawTrack();
  drawReads();
  renderSite();
  if ($("sort").value === "allele") $("sort").onchange();
};
function renderSite() {
  if (!model || site === null) {
    $("sitePanel").hidden = true;
    $("siteLabel").textContent = "No site selected";
    return;
  }
  $("siteLabel").textContent = `${model.region.chr}:${site + 1}`;
  $("sitePanel").hidden = false;
  $("siteTitle").textContent =
    `${model.region.chr}:${site + 1} · REF ${model.reference[site - model.region.start]}`;
  const root = $("siteStats"),
    s = siteCounts(
      model,
      model.fragments.map((f) => f.id),
      site,
    );
  root.replaceChildren(
    ...[...s.counts].map(([b, n]) =>
      el(
        "p",
        `${b}: ${n}/${s.total} callable pairs (${pct(n, s.total)})`,
        "detail-line",
      ),
    ),
    el(
      "p",
      `${s.unknown} uncovered/filtered · ${s.conflicting} contradictory pairs`,
      "muted",
    ),
  );
}

function renderGroups() {
  const root = $("groupList");
  root.replaceChildren();
  if (!model) return;
  const total = model.activeIDs(model.fragments.map((f) => f.id)).length;
  for (const g of model.groups) {
    const members = model.activeIDs(g.members).length,
      button = el("button", undefined, "group-chip");
    button.dataset.group = String(g.id);
    button.append(
      el("span", "●", "swatch"),
      el("strong", g.name),
      el("span", `${members} pairs · ${pct(members, total)}`),
    );
    button.firstChild.style.color = color(g.id);
    button.classList.toggle(
      "selected",
      inspected?.kind === "group" && inspected.id === g.id,
    );
    button.onclick = () => {
      inspect({ kind: "group", id: g.id });
      const i = rows.findIndex(
        (r) => r.kind === "group" && r.group.id === g.id,
      );
      $("readScroll").scrollTop = Math.max(0, packed.byItem.get(i)?.y || 0);
    };
    root.append(button);
  }
}
function inspect(value) {
  inspected = value;
  for (const button of document.querySelectorAll(".group-chip"))
    button.classList.toggle(
      "selected",
      value?.kind === "group" && +button.dataset.group === value.id,
    );
  const root = $("inspector");
  root.replaceChildren();
  if (!model || !value) {
    $("inspectorTitle").textContent = "Select a pair or haplotype";
    root.append(
      el(
        "p",
        "Hover reads for MAPQ. Select pairs to construct a haplotype; groups expose their fractions and variants here.",
        "muted",
      ),
    );
    return;
  }
  if (value.kind === "group") {
    const g = model.groups.find((g) => g.id === value.id);
    if (!g) {
      inspect(null);
      return;
    }
    const stats = groupStatistics(model, g);
    $("inspectorTitle").textContent = g.name;
    const rename = el("div", undefined, "rename"),
      input = el("input"),
      button = el("button", "Rename");
    input.value = g.name;
    input.maxLength = 80;
    input.setAttribute("aria-label", "Rename haplotype");
    button.onclick = () =>
      run(() => model.rename(g.id, input.value), "Haplotype renamed.");
    rename.append(input, button);
    root.append(
      rename,
      el(
        "p",
        `${stats.members} included / ${stats.assigned} assigned pairs · ${stats.hidden} hidden`,
        "detail-line",
      ),
      el(
        "p",
        `Group pair share: ${stats.members}/${stats.total} (${pct(stats.members, stats.total)})`,
        "detail-line",
      ),
      el(
        "p",
        `${stats.summary.unknown} included pairs lack callable sequence at these settings.`,
        "muted",
      ),
    );
    if (stats.summary.conflicts.length)
      root.append(
        el(
          "p",
          `${stats.summary.conflicts.length} conflicting callable sites. Expand and ungroup contrary pairs, or revise evidence settings. No majority override is applied.`,
          "note warning",
        ),
      );
    root.append(
      el(
        "p",
        "Site VAF across all included callable pairs; within-group ALT fraction shown separately. Shares and VAFs are not cellular clone fractions.",
        "note",
      ),
    );
    if (!stats.variants.length)
      root.append(
        el(
          "p",
          "No callable difference from reference in this haplotype.",
          "muted",
        ),
      );
    const fill = (variants, target) => {
      for (const v of variants) {
        const item = el("div", undefined, "variant-stat");
        item.append(
          el("strong", `${v.pos + 1} ${v.ref || "∅"} → ${v.alt || "∅"}`),
          el(
            "div",
            `Site VAF ${v.altCount}/${v.total} (${pct(v.altCount, v.total)})`,
          ),
          el(
            "div",
            `Within group ${v.localAlt}/${v.localTotal} (${pct(v.localAlt, v.localTotal)})`,
          ),
        );
        item.onclick = () => {
          site = v.pos;
          const span = Math.min(100, model.region.end - model.region.start);
          setView(v.pos - span / 2, span);
          renderSite();
        };
        target.append(item);
      }
    };
    fill(stats.variants.slice(0, 12), root);
    if (stats.variants.length > 12) {
      const details = el("details");
      details.append(
        el("summary", `${stats.variants.length - 12} more variant positions`),
      );
      fill(stats.variants.slice(12), details);
      root.append(details);
    }
    root.append(
      el(
        "p",
        `${stats.summary.clips.length} distinct unaligned soft-clip sequences remain available on expansion.`,
        "muted",
      ),
    );
    const expand = el(
      "button",
      model.expanded.has(g.id) ? "Collapse read pairs" : "Expand read pairs",
    );
    expand.onclick = () => {
      if (model.expanded.has(g.id)) model.expanded.delete(g.id);
      else model.expanded.add(g.id);
      rows = model.rows($("sort").value, site, $("showHidden").checked);
      sizeRows();
      drawReads();
      renderAccessible();
      inspect(value);
    };
    const dissolve = el("button", "Dissolve haplotype");
    dissolve.onclick = () =>
      run(() => {
        model.dissolve(g.id);
        inspected = null;
      }, "Haplotype dissolved; every read pair retained.");
    root.append(expand, document.createTextNode(" "), dissolve);
  } else {
    const f = model.fragments[value.id];
    if (!f) return;
    $("inspectorTitle").textContent = f.name;
    for (const rid of f.readIds) {
      const r = model.readById.get(rid);
      root.append(
        el(
          "p",
          `${r.flags & 128 ? "Read 2" : r.flags & 64 ? "Read 1" : "Read"} · MAPQ ${r.mq === 255 ? "unknown" : r.mq} · ${r.flags & 16 ? "reverse" : "forward"}`,
          "detail-line",
        ),
        el(
          "p",
          `${r.chr || model.region.chr}:${r.start + 1}–${r.end} · ${r.cigar || "no CIGAR"} · flags ${r.flags}`,
          "detail-line",
        ),
      );
      if (r.displayIssue) root.append(el("p", r.displayIssue, "note warning"));
      const details = el("details");
      details.append(
        el("summary", "Read sequence"),
        el("div", hasSequence(r) ? r.seq : "Sequence unavailable", "sequence"),
      );
      root.append(details);
      if (r.sa) root.append(el("p", "SA: " + r.sa, "detail-line"));
    }
    root.append(
      el(
        "p",
        "Selection is the whole physical pair, including any available supplementary records.",
        "note",
      ),
    );
  }
}
function drawTree() {
  const canvas = $("tree"),
    height = Math.max(85, (treeData?.groups?.length || 1) * 30 + 30);
  canvas.style.height = height + "px";
  const { c, w } = resize(canvas, height);
  c.fillStyle = "#fff";
  c.fillRect(0, 0, w, height);
  treeHits = [];
  $("cutLabel").textContent = $("treeCut").value + "% SNP difference";
  if (!model || !treeData) return;
  if (treeData.error) {
    $("treeSummary").textContent = treeData.error;
    return;
  }
  const t = treeData;
  $("treeSummary").textContent =
    `${t.groups.length} manual haplotypes → ${t.clusters.length} similarity clusters · ${t.sites.length} variable SNP sites${t.roots.length > 1 ? " · disconnected evidence stays separate" : ""}.`;
  if (!t.groups.length) return;
  const leafX = Math.max(70, w * 0.4),
    scale = leafX - 8,
    pos = new Map();
  let serial = 0;
  const layout = (n) => {
    const y = n.left
      ? (layout(n.left) + layout(n.right)) / 2
      : 25 + serial++ * 30;
    pos.set(n.id, y);
    return y;
  };
  t.roots.forEach(layout);
  const xx = (n) => leafX - n.height * scale;
  c.strokeStyle = "#b1c6d3";
  c.lineWidth = 1.3;
  function paint(n) {
    if (n.left) {
      for (const child of [n.left, n.right]) {
        c.beginPath();
        c.moveTo(xx(n), pos.get(child.id));
        c.lineTo(xx(child), pos.get(child.id));
        c.stroke();
        paint(child);
      }
      c.beginPath();
      c.moveTo(xx(n), pos.get(n.left.id));
      c.lineTo(xx(n), pos.get(n.right.id));
      c.stroke();
    } else {
      const g = model.groups.find((g) => g.id === n.groupID),
        y = pos.get(n.id);
      c.fillStyle = color(g.id);
      c.beginPath();
      c.arc(leafX, y, 3, 0, Math.PI * 2);
      c.fill();
      c.fillStyle = "#456c85";
      c.font = "10px -apple-system,sans-serif";
      c.fillText(g.name.slice(0, 18), leafX + 10, y + 3);
      treeHits.push({ id: g.id, y });
    }
  }
  t.roots.forEach(paint);
  c.setLineDash([3, 3]);
  c.strokeStyle = "#d3af6d";
  const cut = leafX - (+$("treeCut").value / 100) * scale;
  c.beginPath();
  c.moveTo(cut, 5);
  c.lineTo(cut, height);
  c.stroke();
  c.setLineDash([]);
  c.fillStyle = "#859baa";
  c.font = "8px sans-serif";
  c.fillText("SNP similarity · not ancestry", 8, height - 3);
}
$("treeCut").oninput = () => {
  if (model) {
    treeData = manualHierarchy(model, +$("treeCut").value / 100);
    drawTree();
  }
};
$("tree").onclick = (e) => {
  if (!model) return;
  const b = $("tree").getBoundingClientRect(),
    hit = treeHits.find((h) => Math.abs(h.y - (e.clientY - b.top)) < 13);
  if (hit) {
    inspect({ kind: "group", id: hit.id });
    const g = model.groups.find((g) => g.id === hit.id);
    for (const id of g.members)
      if (!model.hidden.has(id) || $("showHidden").checked)
        model.selection.add(id);
    drawReads();
    buttons();
    renderAccessible();
  }
};
$("tree").onmousemove = (e) => {
  if (!model || !hoverDetails) {
    $("tree").title = "";
    return;
  }
  const b = $("tree").getBoundingClientRect(),
    hit = treeHits.find((h) => Math.abs(h.y - (e.clientY - b.top)) < 13);
  if (hit) {
    const g = model.groups.find((g) => g.id === hit.id),
      s = groupStatistics(model, g);
    $("tree").title =
      `${g.name}: ${s.members}/${s.total} pair share (${pct(s.members, s.total)}). ${s.variants
        .slice(0, 6)
        .map(
          (v) =>
            `${v.pos + 1} ${v.alt}: site VAF ${v.altCount}/${v.total} (${pct(v.altCount, v.total)})`,
        )
        .join("; ")}`;
  } else $("tree").title = "";
};
function renderAccessible(page = 0) {
  const root = $("accessibleRows");
  if (!model) {
    root.replaceChildren();
    return;
  }
  const order = packed.order,
    start = page * 50,
    end = Math.min(start + 50, order.length);
  root.replaceChildren(
    el(
      "span",
      order.length
        ? `${start + 1}–${end} of ${order.length} displayed items`
        : "No displayed items",
    ),
  );
  if (start > 0) {
    const b = el("button", "Previous 50");
    b.onclick = () => renderAccessible(page - 1);
    root.append(b);
  }
  if (end < order.length) {
    const b = el("button", "Next 50");
    b.onclick = () => renderAccessible(page + 1);
    root.append(b);
  }
  for (const i of order.slice(start, end)) {
    const row = rows[i],
      p = packed.byItem.get(i),
      b = el(
        "button",
        `${row.kind === "group" ? row.group.name : row.fragment.name}${row.members.some((id) => model.selection.has(id)) ? " ✓" : ""}`,
      );
    b.dataset.item = i;
    b.dataset.x = Math.max(
      leftMargin() + 5,
      Math.min(pile.clientWidth - 18, (p.x1 + p.x2) / 2),
    );
    b.dataset.y = p.y + p.height / 2;
    b.dataset.key = row.key;
    b.setAttribute(
      "aria-pressed",
      String(
        row.members.length > 0 &&
          row.members.every((id) => model.selection.has(id)),
      ),
    );
    b.onclick = (e) => chooseRow(i, e.shiftKey);
    root.append(b);
  }
}

document.onkeydown = (e) => {
  if (/INPUT|SELECT|TEXTAREA/.test(e.target.tagName)) return;
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z") {
    e.preventDefault();
    (e.shiftKey ? $("redo") : $("undo")).click();
  } else if (e.key === "Escape") {
    closeReadMenu();
    $("clearSelection").click();
  } else if (e.shiftKey && e.key.toLowerCase() === "m") {
    $("merge").click();
    e.preventDefault();
  } else if (e.shiftKey && e.key.toLowerCase() === "h") {
    $("hide").click();
    e.preventDefault();
  } else if (e.target === pile && e.key === "ArrowLeft") {
    $("left").click();
    e.preventDefault();
  } else if (e.target === pile && e.key === "ArrowRight") {
    $("right").click();
    e.preventDefault();
  }
};
new ResizeObserver(() => {
  sizeRows();
  drawTrack();
  drawReads();
  drawTree();
}).observe($("readScroll"));
syncChrome();
startWorker();
drawTrack();
drawReads();
buttons();
