import { rawGeometry } from "./alignment-geometry.js";
import { setupHierarchy, groupColor } from "./tree-view.js";
import { consensusGeometry } from "./consensus.js";
import hg19 from "../references/hg19.js";
import hg38 from "../references/hg38.js";
import { ungzip } from "pako";
let catalog;
let catalogBuild;
function genes() {
  const b = $("assembly").value;
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
import { layoutFamilies } from "./families.js";
let expandedFamilies = new Set(),
  layoutCache = null,
  scrollOffset = 0;
const ROW_HEIGHT = 24,
  VIEW_HEIGHT = 480;
const $ = (id) => document.getElementById(id),
  el = (tag, text, cls) => {
    const e = document.createElement(tag);
    if (text !== undefined) e.textContent = text;
    if (cls) e.className = cls;
    return e;
  };
let worker,
  url,
  opened = false,
  data = null,
  view = null,
  rows = 0,
  selected = null,
  job = 0,
  timer,
  contigs = [],
  loci = [],
  bed = [];
const colors = {
  A: "#24865e",
  C: "#3075bc",
  G: "#aa8423",
  T: "#cb575b",
  N: "#b4bdb6",
};
const hierarchyView = setupHierarchy({
  getData: () => data,
  onSelect: (blockID, groupID) => {
    const group = data?.result.overview.find(
      (f) => f.blockID === blockID && f.haplotypeID === groupID,
    );
    if (!group) return;
    $("detail").value = "overview";
    if (
      group.scope &&
      (group.scope.end <= view.start || group.scope.start >= view.end)
    ) {
      const b = data.result.hierarchy.find((b) => b.blockID === blockID),
        span = view.end - view.start,
        center = (b.start + b.end) / 2;
      view.start = Math.max(
        data.region.start,
        Math.min(data.region.end - span, Math.floor(center - span / 2)),
      );
      view.end = view.start + span;
    }
    layoutCache = null;
    const layout = visibleReads(),
      index = layout.findIndex((row) => row.family.id === group.id);
    if (index >= 0) {
      $("readScroll").scrollTop = index * ROW_HEIGHT;
      scrollOffset = $("readScroll").scrollTop;
      draw();
    }
  },
});
function status(text, error = false) {
  $("status").textContent = text;
  $("status").className = error ? "error" : "";
}
function options() {
  const values = {
    mapq: +$("mapq").value,
    baseq: +$("baseq").value,
    support: +$("support").value,
    duplicates: $("duplicates").checked,
    sample: $("sample").value,
  };
  if (
    ![values.mapq, values.baseq, values.support].every(Number.isFinite) ||
    values.mapq < 0 ||
    values.mapq > 254 ||
    values.baseq < 0 ||
    values.baseq > 93 ||
    values.support < 2
  )
    throw Error("Invalid quality or support threshold");
  return values;
}
function busy(on) {
  $("cancel").disabled = !on;
  $("run").disabled = on || !opened;
  $("open").disabled = on;
  $("demo").disabled = on;
  $("demo2").disabled = on;
  clearTimeout(timer);
  if (on)
    timer = setTimeout(() => {
      stop();
      status(
        "Analysis stopped after 60 seconds. Narrow the region before retrying.",
        true,
      );
    }, 60000);
}
function startWorker() {
  worker?.terminate();
  if (url) URL.revokeObjectURL(url);
  url = URL.createObjectURL(
    new Blob([WORKER_SOURCE], { type: "text/javascript" }),
  );
  worker = new Worker(url);
  worker.onerror = () => {
    busy(false);
    status("Analysis worker failed. Clear session and reopen files.", true);
  };
  worker.onmessage = ({ data: m }) => {
    if (m.type === "progress") status(m.message);
    if (m.type === "error") {
      busy(false);
      status(m.message, true);
    }
    if (m.type === "opened") {
      opened = true;
      contigs = m.contigs;
      $("sample").replaceChildren();
      for (const s of m.samples.length ? m.samples : [""])
        $("sample").append(
          new Option(
            s ? "Sample " + ($("sample").options.length + 1) : "Unlabelled",
            s,
          ),
        );
      busy(false);
      status(
        "Files open. Reference identity is user-confirmed; contig lengths will be checked. Choose a region.",
      );
    }
    if (m.type === "result") {
      if (m.request !== job) return;
      busy(false);
      expandedFamilies.clear();
      data = m;
      view = { start: m.region.start, end: m.region.end };
      rows = 0;
      invalidateLayout();
      selected = null;
      render();
      status(
        `${m.result.fragmentCount.toLocaleString()} eligible fragments analyzed. Complete within this interval; no phase inferred outside it. Original alignment records: ${m.result.reads.length.toLocaleString()}.`,
      );
    }
    if (m.type === "event") {
      busy(false);
      $("alignment").textContent =
        `Genomic ref  ${m.path.alignment.a}\nCandidate    ${m.path.alignment.b}\nGenomic ref  ${m.path.alignment.a}`;
      $("eventDetails").textContent =
        `Local anchors ${m.path.start + 1}–${m.path.end}. Dashes in the reference are display columns without genomic positions. Read alignments below are exploratory global affine-gap alignments within bounded windows; original alignments remain available.`;
      $("realigned").textContent = m.realigned
        .map(
          (r) =>
            `Record ${r.id + 1} · score ${r.score}\nCandidate ${r.a}\nRead      ${r.b}`,
        )
        .join("\n\n");
    }
  };
}
function stop() {
  worker?.terminate();
  worker = null;
  opened = false;
  busy(false);
  status(
    "Analysis cancelled. Reopen local files to start a new worker. Previous completed evidence remains visible.",
  );
}
function clear() {
  hierarchyView.clear();
  finishDrag(true);
  suppressFamilyClick = false;
  worker?.terminate();
  worker = null;
  if (url) URL.revokeObjectURL(url);
  url = null;
  clearTimeout(timer);
  opened = false;
  expandedFamilies.clear();
  painted = [];
  layoutCache = null;
  scrollOffset = 0;
  indexedData = null;
  indexedReads = null;
  geometryCache.clear();
  $("readSpacer").style.height = "0px";
  $("readScroll").scrollTop = 0;
  data = null;
  view = null;
  contigs = [];
  loci = [];
  bed = [];
  catalog = null;
  catalogBuild = null;
  job++;
  for (const id of [
    "bam",
    "index",
    "reference",
    "fai",
    "vcf",
    "bed",
    "gene",
    "region",
  ])
    $(id).value = "";
  $("confirm").checked = false;
  $("detail").value = "overview";
  $("collapse").checked = true;
  for (const id of ["locus", "transcript", "sample", "v1", "v2"])
    $(id).replaceChildren();
  for (const id of [
    "blocks",
    "variants",
    "events",
    "structural",
    "readTable",
    "readDetail",
    "alignment",
    "realigned",
    "eventDetails",
    "metrics",
    "pairCounts",
  ])
    $(id).replaceChildren();
  $("pileup").width = 1;
  $("pileup").height = 1;
  delete $("pileup").dataset.viewStart;
  delete $("pileup").dataset.viewEnd;
  $("rowCount").textContent = "";
  $("pileup").removeAttribute("title");
  $("workspace").hidden = true;
  $("expanded").hidden = true;
  $("empty").hidden = false;
  $("regionTitle").textContent = "Your local window into phase";
  busy(false);
  status("Session cleared. Files, analysis, and worker memory released.");
}
$("clear").onclick = clear;
$("cancel").onclick = stop;
window.addEventListener("pagehide", clear);
window.addEventListener("pageshow", (e) => {
  if (e.persisted) clear();
});
$("open").onclick = async () => {
  try {
    for (const id of ["bam", "index", "reference"])
      if (!$(id).files[0])
        throw Error("Select BAM, index, and reference files");
    if (!$("confirm").checked)
      throw Error("Confirm the exact alignment reference before opening");
    expandedFamilies.clear();
    layoutCache = null;
    scrollOffset = 0;
    indexedData = null;
    indexedReads = null;
    geometryCache.clear();
    data = null;
    $("workspace").hidden = true;
    $("empty").hidden = false;
    opened = false;
    startWorker();
    busy(true);
    status("Opening selected local files…");
    worker.postMessage({
      type: "open",
      files: Object.fromEntries(
        ["bam", "index", "reference", "fai", "vcf"].map((id) => [
          id,
          $(id).files[0],
        ]),
      ),
    });
  } catch (e) {
    busy(false);
    status(e.message, true);
  }
};
function region() {
  const m = $("region")
    .value.replaceAll(",", "")
    .trim()
    .match(/^([^:]+):(\d+)-(\d+)$/);
  if (!m)
    throw Error("Enter contig:start-end using 1-based inclusive positions");
  let chr = m[1];
  if (contigs.length && !contigs.some((c) => c.refName === chr)) {
    const aliases = contigs.filter(
      (c) => c.refName.replace(/^chr/, "") === chr.replace(/^chr/, ""),
    );
    if (aliases.length === 1) chr = aliases[0].refName;
  }
  return { chr, start: +m[2] - 1, end: +m[3] };
}
$("run").onclick = async () => {
  try {
    if (!worker || !opened) throw Error("Open files first");
    const r = region(),
      opts = options();
    bed = [];
    const bf = $("bed").files[0];
    if (bf) {
      if (bf.size > 8e6) throw Error("BED exceeds 8 MiB");
      for (const line of (await bf.text()).split(/\r?\n/)) {
        if (!line || line.startsWith("#") || line.startsWith("track")) continue;
        const [c, s, e] = line.split("\t");
        if (c.replace(/^chr/, "") === r.chr.replace(/^chr/, "")) {
          if (
            !Number.isSafeInteger(+s) ||
            !Number.isSafeInteger(+e) ||
            +s < 0 ||
            +e <= +s
          )
            throw Error("Invalid BED interval");
          bed.push([+s, +e]);
        }
      }
    }
    busy(true);
    job++;
    worker.postMessage({
      type: "analyze",
      region: r,
      opts,
      request: job,
      mask: $("bedRole").value === "primer" ? bed : [],
    });
  } catch (e) {
    busy(false);
    status(e.message, true);
  }
};
function demo() {
  try {
    clear();
    startWorker();
    busy(true);
    status("Analyzing synthetic SNVs, insertion/duplication, and deletion…");
    worker.postMessage({ type: "demo", opts: options(), request: ++job });
  } catch (e) {
    status(e.message, true);
  }
}
$("demo").onclick = demo;
$("demo2").onclick = demo;
$("lookup").onclick = () => {
  const name = $("gene").value.trim().toUpperCase();
  loci = genes().filter((g) => g[0].toUpperCase() === name);
  $("locus").replaceChildren(new Option("Select a locus", ""));
  loci.forEach((g, i) =>
    $("locus").append(new Option(`${g[0]} · ${g[1]}:${g[3] + 1}–${g[4]}`, i)),
  );
  if (!loci.length)
    status(
      "No exact gene symbol in the offline catalog. Try an official symbol or enter genomic coordinates.",
      true,
    );
  else if (loci.length === 1) {
    $("locus").value = "0";
    $("locus").onchange();
  }
};
$("locus").onchange = () => {
  const g = loci[$("locus").value];
  if (!g) return;
  $("region").value = `${g[1]}:${g[3] + 1}-${g[4]}`;
  $("transcript").replaceChildren(new Option("Whole gene", ""));
  g[6].forEach((t, i) => $("transcript").append(new Option(t[0], i)));
};
$("transcript").onchange = () => {
  const g = loci[$("locus").value];
  if (!g) return;
  const t = $("transcript").value === "" ? null : g[6][+$("transcript").value];
  $("region").value = `${g[1]}:${(t ? t[1] : g[3]) + 1}-${t ? t[2] : g[4]}`;
};
function visibleReads() {
  if (!data) return [];
  if (!layoutCache)
    layoutCache = layoutFamilies(displayGroups(), data.result.reads, {
      start: view.start,
      end: view.end,
      layer: $("layer").value,
      sort: $("sort").value,
      collapsed: $("collapse").checked,
      expanded: expandedFamilies,
    });
  return layoutCache;
}
function displayGroups() {
  return $("detail").value === "overview"
    ? data.result.overview
    : data.result.families;
}
function invalidateLayout() {
  layoutCache = null;
  scrollOffset = 0;
  $("readScroll").scrollTop = 0;
}
function toggleFamily(id) {
  if (expandedFamilies.has(id)) expandedFamilies.delete(id);
  else expandedFamilies.add(id);
  layoutCache = null;
  draw();
}
function render() {
  const r = data.result;
  $("empty").hidden = true;
  $("workspace").hidden = false;
  $("regionTitle").textContent =
    `${data.region.chr}:${data.region.start + 1}–${data.region.end}`;
  $("metrics").replaceChildren();
  for (const [v, label] of [
    [r.fragmentCount, "Eligible fragments"],
    [r.variants.length, "Candidate alleles"],
    [r.blocks.length, "Local phase blocks"],
    [r.excluded, "Excluded records"],
  ]) {
    const d = el("div", undefined, "metric");
    d.append(el("b", v.toLocaleString()), el("span", label));
    $("metrics").append(d);
  }
  $("blocks").replaceChildren();
  for (const b of r.blocks) {
    const div = el("div", undefined, "block");
    div.append(el("b", `Block ${b.id} · ${b.start + 1}–${b.end}`));
    for (const g of b.groups) {
      const d = el("div", undefined, "group");
      d.append(
        el("b", `H${g.id} · ${g.count} assigned fragments`),
        el(
          "div",
          g.pattern
            .map((a, i) => `${r.variants[b.indices[i]].pos + 1}:${a}`)
            .join(" · "),
        ),
        el("small", `${g.compatible} compatible but ambiguous`),
      );
      div.append(d);
    }
    div.append(
      el(
        "p",
        `${b.ambiguous} ambiguous · ${b.conflicting} unsupported/conflicting. ${b.indices.length === 1 ? "Single site: no inter-variant phase." : "Links require repeated fragment co-observation."}`,
      ),
    );
    $("blocks").append(div);
  }
  if (!r.blocks.length)
    $("blocks").append(
      el("p", "No candidate alleles meet the current filters."),
    );
  $("variants").replaceChildren();
  for (const v of r.variants) {
    const row = el("tr");
    for (const x of [
      v.pos + 1,
      `${v.ref || "∅"} → ${v.alt || "∅"}`,
      v.support,
      v.overlap
        ? v.complexOverlap
          ? "Overlapping allele · inspect locally"
          : "Multiallelic SNV"
        : v.supplied
          ? "Supplied / observed"
          : "Discovered",
    ])
      row.append(el("td", x));
    row.onclick = () => {
      view = {
        start: Math.max(data.region.start, v.pos - 70),
        end: Math.min(data.region.end, v.pos + 100),
      };
      rows = 0;
      invalidateLayout();
      draw();
    };
    $("variants").append(row);
  }
  for (const id of ["v1", "v2"]) {
    $(id).replaceChildren();
    r.variants.forEach((v, i) =>
      $(id).append(
        new Option(`${v.pos + 1} ${v.ref || "∅"}>${v.alt || "∅"}`, i),
      ),
    );
  }
  if (r.variants.length > 1) $("v2").value = "1";
  pair();
  $("events").replaceChildren();
  for (const e of r.events) {
    const b = el("button", `${e.type} · ${e.pos + 1} · ${e.support} fragments`);
    b.onclick = () => showEvent(e);
    $("events").append(b);
  }
  for (const c of r.clips)
    $("events").append(
      el(
        "p",
        `Unresolved ${c.side} soft clip at ${c.pos + 1}: ${c.seq.length} bases · ${c.support} fragments. Shared clip sequence alone does not identify an event.`,
        "muted",
      ),
    );
  if (!r.events.length)
    $("events").append(
      el("p", "No supported insertion or deletion candidates."),
    );
  renderStructural(r.structural);
  hierarchyView.render();
  draw();
}
function renderStructural(structural) {
  const panel = $("structural");
  panel.replaceChildren();
  for (const warning of structural?.warnings || [])
    panel.append(el("p", warning, "muted"));
  if (!structural?.events.length) {
    panel.append(
      el(
        "p",
        "No split-alignment or discordant-mate signals found in the available local records. Missing SA/mate metadata can limit detection.",
        "muted",
      ),
    );
    return;
  }
  for (const event of structural.events) {
    const card = el("details", undefined, "block");
    card.append(
      el(
        "summary",
        `${event.type} · ${event.fragments} fragment${event.fragments === 1 ? "" : "s"} · ${event.partner.chr}:${event.partnerRange[0] + 1}`,
      ),
    );
    card.append(
      el(
        "p",
        `Local ${event.local.chr}:${event.localRange[0] + 1}–${event.localRange[1] + 1} (${event.local.strand}) → partner ${event.partner.chr}:${event.partnerRange[0] + 1}–${event.partnerRange[1] + 1} (${event.partner.strand}). Anchors clustered within ${structural.clusterTolerance} bp; mate-only coordinates are alignment starts, not exact breakpoints.`,
      ),
    );
    card.append(
      el(
        "p",
        `${event.eligibleFragments} fragments pass mapping-quality filters at both ends; ${event.localQualityFragments} pass local filters; ${event.unknownPartnerQuality} have unknown partner MAPQ; ${event.duplicateFragments} are duplicate-marked. Sources: ${event.sources.join(", ")}. Record IDs: ${event.records.map((id) => id + 1).join(", ")}.`,
      ),
    );
    for (const note of event.notes) card.append(el("p", note, "muted"));
    const button = el(
        "button",
        `Look up partner genes in ${$("assembly").value}`,
      ),
      description = el("p", undefined, "muted");
    button.onclick = () => {
      const matches = genes()
        .filter(
          (g) =>
            g[1].replace(/^chr/, "") ===
              event.partner.chr.replace(/^chr/, "") &&
            g[3] < event.partner.end &&
            g[4] > event.partner.start,
        )
        .map((g) => g[0]);
      description.textContent = matches.length
        ? `Partner alignment overlaps: ${[...new Set(matches)].join(", ")} (${$("assembly").value} catalog). Overlap alone does not establish a gene fusion.`
        : "No overlapping gene in the selected offline assembly catalog. Verify the catalog matches your BAM reference.";
    };
    card.append(button, description);
    panel.append(card);
  }
}
function pair() {
  if (!data) return;
  const a = +$("v1").value,
    b = +$("v2").value;
  const counts = new Map();
  let unknown = 0;
  for (const v of data.result.vectors) {
    if (
      v[a] === "?" ||
      v[b] === "?" ||
      v[a] === undefined ||
      v[b] === undefined
    ) {
      unknown++;
      continue;
    }
    const key = v[a] + " / " + v[b];
    counts.set(key, (counts.get(key) || 0) + 1);
  }
  $("pairCounts").replaceChildren();
  for (const [k, n] of counts)
    $("pairCounts").append(el("p", `${k}: ${n} fragments`));
  $("pairCounts").append(
    el(
      "p",
      `${unknown} fragments not informative at both sites. Absence of joint alternate observations does not establish trans.`,
      "muted",
    ),
  );
}
$("v1").onchange = pair;
$("v2").onchange = pair;
function showEvent(e) {
  $("alignment").textContent = "Aligning every overlapping record…";
  $("realigned").textContent = "";
  $("eventDetails").textContent = "";
  $("expanded").hidden = false;
  $("eventTitle").textContent =
    `${e.type} at ${e.pos + 1}${e.source ? ` · source ${e.source[0] + 1}–${e.source[1]}` : ""}`;
  if (e.alt.length > 1000 || e.ref.length > 5000) {
    $("alignment").textContent =
      "Event exceeds local reconstruction limits. Inspect the original pileup.";
    return;
  }
  if (!worker) {
    status("Reopen files before requesting local alignment.", true);
    return;
  }
  busy(true);
  worker.postMessage({
    type: "event",
    event: e,
    reference: data.reference,
    origin: data.region.start,
    reads: data.result.reads,
  });
}
let painted = [];
function draw() {
  if (!data) return;
  $("pileup").dataset.viewStart = String(view.start);
  $("pileup").dataset.viewEnd = String(view.end);
  const canvas = $("pileup"),
    width = Math.max(400, $("canvasWrap").clientWidth),
    all = visibleReads();
  const totalHeight = all.length * ROW_HEIGHT;
  const viewportHeight = Math.min(VIEW_HEIGHT, Math.max(120, totalHeight));
  $("readSpacer").style.height = Math.max(viewportHeight, totalHeight) + "px";
  scrollOffset = Math.min(
    scrollOffset,
    Math.max(0, totalHeight - viewportHeight),
  );
  rows = Math.floor(scrollOffset / ROW_HEIGHT);
  painted = all.slice(rows, rows + Math.ceil(viewportHeight / ROW_HEIGHT) + 1);
  const height = 110 + viewportHeight,
    dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr;
  canvas.height = height * dpr;
  canvas.style.height = height + "px";
  const c = canvas.getContext("2d");
  c.scale(dpr, dpr);
  c.fillStyle = "#fff";
  c.fillRect(0, 0, width, height);
  const left = 115,
    scale = (width - left - 15) / (view.end - view.start),
    x = (p) => left + (p - view.start) * scale;
  c.font = "10px monospace";
  c.fillStyle = "#648073";
  c.fillText("REFERENCE", 10, 34);
  c.fillText("COVERAGE", 10, 70);
  const span = view.end - view.start,
    step = Math.max(1, 10 ** Math.floor(Math.log10(span / 6)));
  for (let p = Math.ceil(view.start / step) * step; p < view.end; p += step) {
    c.fillText((p + 1).toLocaleString(), x(p), 14);
    c.strokeStyle = "#edf2ee";
    c.beginPath();
    c.moveTo(x(p), 18);
    c.lineTo(x(p), height);
    c.stroke();
  }
  const max = data.result.coverage.reduce((a, b) => Math.max(a, b), 1);
  for (let p = Math.floor(view.start); p < view.end; p++) {
    const b = data.reference[p - data.region.start];
    c.fillStyle = colors[b] || "#888";
    if (scale > 8) c.fillText(b || "N", x(p), 34);
    const v = data.result.coverage[p - data.region.start] || 0;
    c.fillStyle = "#a5c5b2";
    c.fillRect(x(p), 85 - (v / max) * 35, Math.max(1, scale), (v / max) * 35);
  }
  for (const [s, e] of bed) {
    c.fillStyle = "#e3b56855";
    c.fillRect(x(s), 40, (e - s) * scale, 45);
  }
  for (const b of data.result.blocks) {
    c.strokeStyle = "#47745f";
    c.beginPath();
    c.moveTo(x(b.start), 94);
    c.lineTo(x(b.end), 94);
    c.stroke();
  }
  c.save();
  c.beginPath();
  c.rect(0, 110, width, VIEW_HEIGHT);
  c.clip();
  painted.forEach((row, i) => {
    const y = 110 + i * ROW_HEIGHT - (scrollOffset % ROW_HEIGHT),
      family = row.family;
    c.fillStyle = row.type === "family" ? "#edf4ef" : "#fff";
    c.fillRect(0, y, width, ROW_HEIGHT - 1);
    const stats = data.result.hierarchy
      ?.find((b) => b.blockID === family.blockID)
      ?.groups.find((g) => g.id === family.haplotypeID);
    if (stats) {
      c.fillStyle = groupColor(family.blockID, family.haplotypeID);
      c.fillRect(0, y, 3, ROW_HEIGHT - 1);
    }
    c.fillStyle = "#365f50";
    const label =
      row.type === "family"
        ? `${expandedFamilies.has(family.id) ? "▾" : "▸"} ${family.kind === "overview" ? family.label : `F${family.id + 1}`} ×${row.count}`
        : `  ${row.member.paired ? (row.member.matePresent ? "Pair" : "Single end") : "Read"} ${row.member.id + 1}`;
    c.fillText(
      label +
        (family.discordantPair ? " · !" : "") +
        (family.structuralLinks?.length ? " · SV" : ""),
      6,
      y + 10,
      left - 12,
    );
    c.fillStyle = "#7d8e84";
    c.fillText(
      stats
        ? `${stats.fragments} assigned · ${(100 * (stats.assignedShare || 0)).toFixed(1)}%`
        : row.group.slice(0, 17),
      6,
      y + 21,
      left - 12,
    );
    c.save();
    c.beginPath();
    c.rect(left, y, width - left - 20, ROW_HEIGHT);

    c.clip();
    if (family.scope) {
      c.beginPath();
      c.rect(
        x(family.scope.start),
        y,
        Math.max(0, (family.scope.end - family.scope.start) * scale),
        ROW_HEIGHT,
      );
      c.clip();
    }
    const visibleIDs =
      row.type === "family"
        ? row.members.flatMap((m) => m.readIds)
        : row.member.readIds;
    const members = row.type === "family" ? row.members : [row.member];
    // Draw each fragment's connector, never a family-wide fabricated bridge.
    c.strokeStyle = "#b6c5bc";
    c.setLineDash([3, 3]);
    const byID = readIndex();
    for (const member of members) {
      const rr = member.readIds.map((id) => byID.get(id)),
        r1 = rr.find((r) => r.flags & 64),
        r2 = rr.find((r) => r.flags & 128);
      if (r1 && r2) {
        c.beginPath();
        c.moveTo(x(Math.min(r1.end, r2.end)), y + 11);
        c.lineTo(x(Math.max(r1.start, r2.start)), y + 11);
        c.stroke();
      }
    }
    c.setLineDash([]);
    const glyph = geometry(row, visibleIDs);
    c.fillStyle = visibleIDs.every((id) => byID.get(id).excluded)
      ? "#dedbd5"
      : "#b9cbbf";
    for (const [s, e] of glyph.segments)
      c.fillRect(x(s), y + 6, Math.max(1, (e - s) * scale), 10);
    for (const m of glyph.mismatches) {
      if (m.pos < view.start || m.pos >= view.end) continue;
      c.fillStyle = colors[m.base] || "#aaa";
      c.fillRect(x(m.pos), y + 6, Math.max(1, scale), 10);
      if (scale > 9) {
        c.fillStyle = "#fff";
        c.fillText(m.base, x(m.pos), y + 15);
      }
    }
    for (const conflict of glyph.conflicts || []) {
      const xp = x(conflict.pos);
      c.strokeStyle = "#b66e18";
      c.beginPath();
      c.moveTo(xp, y + 5);
      c.lineTo(xp + Math.max(4, scale), y + 17);
      c.moveTo(xp, y + 17);
      c.lineTo(xp + Math.max(4, scale), y + 5);
      c.stroke();
    }
    for (const e of glyph.events) {
      if (e.op === "I") {
        c.fillStyle = "#9c70b3";
        c.fillRect(x(e.rawPos ?? e.pos) - 2, y + 3, 4, 16);
      } else {
        c.strokeStyle = e.op === "D" ? "#7b817c" : "#bcc9c0";
        c.beginPath();
        c.moveTo(x(e.rawPos ?? e.pos), y + 11);
        c.lineTo(x((e.rawPos ?? e.pos) + e.length), y + 11);
        c.stroke();
      }
    }
    for (const clip of glyph.clips) {
      c.fillStyle = "#b3a0c2";
      c.fillRect(x(clip.pos) - 2, y + 6, 4, 10);
      if ($("clips").checked && scale > 8)
        c.fillText(clip.seq, x(clip.pos), y + 15);
    }
    c.restore();
  });
  c.restore();
  const families = new Set(all.map((row) => row.family.id));
  const fragmentCount = data.result.families.reduce((n, f) => n + f.count, 0);
  $("rowCount").textContent =
    `${families.size} ${$("detail").value === "overview" ? "group rows" : "families"} in view · rows ${all.length ? rows + 1 : 0}–${Math.min(rows + Math.ceil(viewportHeight / ROW_HEIGHT), all.length)} of ${all.length} · ${fragmentCount} total fragments / ${data.result.reads.length} records retained · no downsampling`;
  $("readTable").replaceChildren();
  for (const row of painted) {
    const tr = el("tr"),
      family = row.family,
      member = row.type === "fragment" ? row.member : null;
    for (const value of [
      row.type === "family"
        ? `${family.kind === "overview" ? family.label : `Family ${family.id + 1}`} · ${row.count} fragments`
        : `${member.paired ? (member.matePresent ? "Pair" : "Single end") : "Read"} ${member.id + 1}`,
      row.start + 1,
      member
        ? member.readIds.map((id) => readIndex().get(id).cigar).join(" / ")
        : `${family.recordCount} records · ${family.duplicates} duplicate-marked fragments`,
      member
        ? member.readIds.map((id) => readIndex().get(id).mq).join(" / ")
        : "Expand to inspect",
      member
        ? member.readIds.map((id) => readIndex().get(id).flags).join(" / ")
        : `${family.excluded} excluded fragments`,
      row.group,
    ])
      tr.append(el("td", value));
    tr.tabIndex = 0;
    tr.setAttribute("role", "button");
    if (row.type === "family")
      tr.setAttribute("aria-expanded", String(expandedFamilies.has(family.id)));
    const detail = () => {
      if (row.type === "family") {
        toggleFamily(family.id);
        return;
      }
      $("readDetail").textContent = member.readIds
        .map((id) => {
          const r = readIndex().get(id);
          return `Read name (selected explicitly): ${r.name}\nStart: ${r.start + 1}\nCIGAR: ${r.cigar}\nFlags: ${r.flags}\nMate: ${r.mateChr || "unknown contig"}:${r.mate >= 0 ? r.mate + 1 : "unavailable"}\nSA: ${r.sa || "not present"}\nSequence: ${r.seq}\nBase qualities: ${r.qual.join(",")}`;
        })
        .join("\n\n");
    };
    tr.onclick = detail;
    tr.onkeydown = (e) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        detail();
      }
    };
    $("readTable").append(tr);
  }
}
let indexedData = null,
  indexedReads = null,
  geometryCache = new Map();
function readIndex() {
  if (indexedData !== data) {
    indexedData = data;
    indexedReads = new Map(data.result.reads.map((r) => [r.id, r]));
    geometryCache.clear();
  }
  return indexedReads;
}
function geometry(row, ids) {
  const key =
    row.type === "family"
      ? `${row.family.kind === "overview" ? "g" : "f"}${row.family.id}:${$("layer").value}`
      : `m${row.member.id}:${$("layer").value}`;
  if (geometryCache.has(key)) return geometryCache.get(key);
  const summaryMembers =
    row.type === "family"
      ? row.family.members
          .map((m) => ({
            ...m,
            readIds: m.readIds.filter(
              (id) =>
                $("layer").value === "all" ||
                ($("layer").value === "excluded"
                  ? readIndex().get(id).excluded
                  : !readIndex().get(id).excluded),
            ),
          }))
          .filter((m) => m.readIds.length)
      : null;
  if (summaryMembers) ids = summaryMembers.flatMap((m) => m.readIds);

  if (row.type === "family" && row.family.kind === "overview") {
    const glyph = consensusGeometry(
      summaryMembers,
      readIndex(),
      data.reference,
      data.region.start,
      data.result.settings,
    );
    geometryCache.set(key, glyph);
    return glyph;
  }
  const g = rawGeometry(
    ids.map((id) => readIndex().get(id)),
    data.reference,
    data.region.start,
  );
  geometryCache.set(key, g);
  return g;
}
let dragState = null,
  dragFrame = 0,
  suppressFamilyClick = false;
function applyDrag() {
  if (!dragState?.moved || !data || !view) return;
  const span = dragState.end - dragState.start;
  const next = Math.max(
    data.region.start,
    Math.min(
      data.region.end - span,
      Math.round(
        dragState.start +
          ((dragState.x - dragState.currentX) * span) / dragState.width,
      ),
    ),
  );
  if (view.start === next) return;
  view.start = next;
  view.end = next + span;
  layoutCache = null;
  draw();
}
function finishDrag(cancelled = false) {
  if (dragFrame) cancelAnimationFrame(dragFrame);
  dragFrame = 0;
  if (!dragState) return;
  if (!cancelled) applyDrag();
  const id = dragState.id;
  suppressFamilyClick = dragState.moved;
  dragState = null;
  const canvas = $("pileup");
  canvas.classList.remove("dragging");
  if (canvas.hasPointerCapture(id)) canvas.releasePointerCapture(id);
}
$("pileup").onpointerdown = (e) => {
  if (!data || !view || e.button !== 0 || !e.isPrimary || dragState) return;
  const canvas = $("pileup"),
    bounds = canvas.getBoundingClientRect();
  suppressFamilyClick = false;
  dragState = {
    id: e.pointerId,
    x: e.clientX,
    currentX: e.clientX,
    start: view.start,
    end: view.end,
    width: Math.max(1, bounds.width - 130),
    moved: false,
  };
  canvas.setPointerCapture(e.pointerId);
  canvas.focus({ preventScroll: true });
};
$("pileup").onpointermove = (e) => {
  if (!dragState || e.pointerId !== dragState.id) return;
  dragState.currentX = e.clientX;
  if (Math.abs(e.clientX - dragState.x) > 4) {
    dragState.moved = true;
    $("pileup").classList.add("dragging");
  }
  if (dragState.moved && !dragFrame)
    dragFrame = requestAnimationFrame(() => {
      dragFrame = 0;
      applyDrag();
    });
};
$("pileup").onpointerup = (e) => {
  if (dragState?.id === e.pointerId) finishDrag();
};
$("pileup").onpointercancel = (e) => {
  if (dragState?.id === e.pointerId) finishDrag(true);
};
$("pileup").onlostpointercapture = (e) => {
  if (dragState?.id === e.pointerId) finishDrag(true);
};
$("pileup").addEventListener("pointermove", (e) => {
  if (dragState || !data) return;
  const y = e.clientY - $("pileup").getBoundingClientRect().top;
  if (y < 110) {
    $("pileup").removeAttribute("title");
    return;
  }
  const row = visibleReads()[Math.floor((y - 110 + scrollOffset) / ROW_HEIGHT)];
  if (!row) return;
  const f = row.family;
  if (f.blockID && f.haplotypeID) {
    hierarchyView.showGroup(f.blockID, f.haplotypeID);
    const b = data.result.hierarchy.find((b) => b.blockID === f.blockID),
      g = b.groups.find((g) => g.id === f.haplotypeID);
    $("pileup").title = hierarchyView.info(b, g);
  } else {
    const text = `${f.label || "Family"}: ${row.count} displayed fragments; ${f.excluded} excluded, ${f.duplicates} duplicate-marked. These records are not assigned to a unique supported haplotype. Expand to inspect every pair.`;
    hierarchyView.showSummary(text);
    $("pileup").title = text;
  }
});
$("pileup").onclick = (e) => {
  if (suppressFamilyClick) {
    suppressFamilyClick = false;
    return;
  }
  const y = e.clientY - $("pileup").getBoundingClientRect().top;
  if (y < 110) return;
  const row = visibleReads()[Math.floor((y - 110 + scrollOffset) / ROW_HEIGHT)];
  if (row?.type === "family") toggleFamily(row.family.id);
};
$("readScroll").onscroll = () => {
  scrollOffset = $("readScroll").scrollTop;
  draw();
};
$("pileup").addEventListener(
  "wheel",
  (e) => {
    e.preventDefault();
    $("readScroll").scrollTop += e.deltaY;
  },
  { passive: false },
);
$("collapse").onchange = () => {
  invalidateLayout();
  draw();
};
$("collapseAll").onclick = () => {
  expandedFamilies.clear();
  $("collapse").checked = true;
  invalidateLayout();
  draw();
};
$("expandAll").onclick = () => {
  for (const f of data ? displayGroups() : []) expandedFamilies.add(f.id);
  invalidateLayout();
  draw();
};

function zoom(f) {
  if (!view) return;
  const center = (view.start + view.end) / 2,
    span = Math.max(
      20,
      Math.min(
        data.region.end - data.region.start,
        (view.end - view.start) * f,
      ),
    );
  view.start = Math.max(
    data.region.start,
    Math.min(data.region.end - span, Math.floor(center - span / 2)),
  );
  view.end = Math.min(data.region.end, view.start + Math.ceil(span));
  rows = 0;
  invalidateLayout();
  draw();
}
function pan(f) {
  if (!view) return;
  const span = view.end - view.start;
  view.start = Math.max(
    data.region.start,
    Math.min(data.region.end - span, Math.round(view.start + span * f)),
  );
  view.end = view.start + span;
  rows = 0;
  invalidateLayout();
  draw();
}
$("zoomIn").onclick = () => zoom(0.5);
$("zoomOut").onclick = () => zoom(2);
$("left").onclick = () => pan(-0.5);
$("right").onclick = () => pan(0.5);
$("fit").onclick = () => {
  if (data) {
    view = { ...data.region };
    rows = 0;
    invalidateLayout();
    draw();
  }
};
for (const id of ["sort", "layer", "clips"])
  $(id).onchange = () => {
    rows = 0;
    invalidateLayout();
    draw();
  };
$("detail").onchange = () => {
  expandedFamilies.clear();
  geometryCache.clear();
  invalidateLayout();
  draw();
};
$("pileup").onkeydown = (e) => {
  if (
    ["ArrowDown", "ArrowUp", "PageDown", "PageUp", "Home", "End"].includes(
      e.key,
    )
  ) {
    e.preventDefault();
    const scroll = $("readScroll");
    if (e.key === "Home") scroll.scrollTop = 0;
    else if (e.key === "End") scroll.scrollTop = scroll.scrollHeight;
    else
      scroll.scrollTop +=
        (e.key === "ArrowUp"
          ? -1
          : e.key === "ArrowDown"
            ? 1
            : e.key === "PageUp"
              ? -20
              : 20) * ROW_HEIGHT;
    return;
  }
  if (e.key === "Enter") {
    const row = visibleReads()[Math.floor(scrollOffset / ROW_HEIGHT)];
    if (row?.type === "family") toggleFamily(row.family.id);
    return;
  }
  if (["ArrowLeft", "ArrowRight", "+", "-"].includes(e.key)) {
    e.preventDefault();
    if (e.key === "ArrowLeft") pan(-0.25);
    else if (e.key === "ArrowRight") pan(0.25);
    else zoom(e.key === "+" ? 0.5 : 2);
  }
};
new ResizeObserver(() => {
  draw();
  hierarchyView.draw();
}).observe($("canvasWrap"));
