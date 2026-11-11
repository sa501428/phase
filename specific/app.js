import { watchVaf } from "../src/live-vaf.js";
import { pairMetrics } from "./phasing.js";
import { helperGraph, findHelperPath } from "./helpers.js";
const $ = (id) => document.getElementById(id);
const node = (tag, text, cls) => {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  if (cls) e.className = cls;
  return e;
};
const fmt = (n) => n.toLocaleString();
const percent = (n) => (n === null ? "—" : (100 * n).toFixed(1) + "%");
let worker,
  workerURL,
  timer,
  busy = false,
  opened = false,
  parsed = null,
  result = null,
  pairMap,
  graph,
  components;
let parsePage = 0,
  pairPage = 0,
  helperPage = 0,
  selected = null,
  usedHelpers = new Set();
let requestedPairs = [],
  totalAll = 0;
function status(message, error = false) {
  $("status").textContent = message;
  $("status").classList.toggle("error", error);
}
function controls() {
  for (const id of ["open", "parse", "run"])
    $(id).disabled =
      busy ||
      (id === "run" &&
        (!opened ||
          !parsed?.variants.length ||
          parsed.entries.some((e) => e.error)));
  for (const id of [
    "bam",
    "index",
    "variants",
    "assembly",
    "sample",
    "mapq",
    "baseq",
    "duplicates",
    "helpers",
    "helperSupport",
    "helperFraction",
    "helperGap",
    "linkAgreement",
  ])
    $(id).disabled = busy || (id === "sample" && !opened);
  $("cancel").hidden = !busy;
  for (const e of $("helperSettings").querySelectorAll("input"))
    e.disabled = busy || !$("helpers").checked;
}
function finish() {
  busy = false;
  clearTimeout(timer);
  controls();
}
function armWatchdog() {
  clearTimeout(timer);
  timer = setTimeout(() => {
    stopWorker();
    opened = false;
    invalidateParsed();
    clearResults();
    finish();
    status(
      "Processing exceeded three minutes without progress and was stopped. Reopen files and reduce the target list or helper search scope. No partial result was reported.",
      true,
    );
  }, 180000);
}
function stopWorker() {
  clearTimeout(timer);
  worker?.terminate();
  worker = null;
  if (workerURL) URL.revokeObjectURL(workerURL);
  workerURL = null;
}
function startWorker() {
  if (worker) return;
  workerURL = URL.createObjectURL(
    new Blob([WORKER_SOURCE], { type: "text/javascript" }),
  );
  worker = new Worker(workerURL);
  worker.onerror = (event) => {
    stopWorker();
    opened = false;
    finish();
    status(event.message || "Worker failed; reopen files.", true);
  };
  worker.onmessage = ({ data: m }) => {
    if (m.type === "progress") {
      status(m.message);
      armWatchdog();
    } else if (m.type === "error") {
      finish();
      status(m.message, true);
    } else if (m.type === "opened") {
      opened = true;
      $("sample").replaceChildren();
      for (const name of m.samples) {
        const o = node("option", name);
        o.value = name;
        $("sample").append(o);
      }
      $("fileStatus").textContent =
        `${fmt(m.contigs)} contigs · ${m.samples.length} sample${m.samples.length === 1 ? "" : "s"}`;
      finish();
      status("Opened indexed BAM. Choose a sample and check the variant list.");
      if ($("variants").value.trim()) checkVariants();
    } else if (m.type === "parsed") {
      parsed = m;
      parsePage = 0;
      renderParsed();
      finish();
      const errors = m.entries.filter((e) => e.error).length;
      status(
        errors
          ? `${errors} variant line${errors === 1 ? "" : "s"} need correction. All errors must be resolved before analysis.`
          : `Checked ${fmt(m.variants.length)} unique variant${m.variants.length === 1 ? "" : "s"}. ${opened ? "Ready to phase." : "Open BAM/BAI to continue."}`,
        !!errors,
      );
    } else if (m.type === "result") {
      result = m.result;
      finish();
      renderResult();
      status(
        `Complete · ${fmt(result.requestedCount)} requested variants · ${fmt(result.fragments)} read pairs / fragments · ${fmt(result.helperCount)} qualifying helper SNPs. ${result.suppressedVariants?.length || 0} variants suppressed below minimum VAF. Results remain in memory only.`,
      );
    } else if (
      m.type === "evidence" &&
      selected &&
      selected[0] === m.a &&
      selected[1] === m.b
    )
      renderEvidence(m);
  };
}
function clearResults() {
  result = null;
  graph = pairMap = components = null;
  selected = null;
  requestedPairs = [];
  usedHelpers.clear();
  $("results").hidden = true;
  for (const id of [
    "pairRows",
    "helperRows",
    "evidenceRows",
    "pairDetail",
    "pathDetail",
    "summary",
    "variantA",
    "variantB",
  ])
    $(id).replaceChildren();
}
function invalidateParsed() {
  parsed = null;
  $("parseSummary").textContent = "Check the variant list before running";
  $("parseResults").replaceChildren();
  $("parseResults").hidden = true;
  $("parsePager").replaceChildren();
  $("parsePager").hidden = true;
}
function reset() {
  stopWorker();
  opened = false;
  busy = false;
  clearResults();
  invalidateParsed();
  for (const id of ["bam", "index", "variants"]) $(id).value = "";
  $("sample").replaceChildren(node("option", "Open BAM first"));
  $("sample").firstChild.value = "";
  $("fileStatus").textContent = "No files opened";
  controls();
  status("Session cleared. No case data is retained by this app.");
}
function pager(id, count, page, size, callback) {
  const target = $(id),
    pages = Math.ceil(count / size);
  target.replaceChildren();
  target.hidden = !count;
  if (!count) return;
  const prev = node("button", "← Previous"),
    next = node("button", "Next →");
  prev.type = next.type = "button";
  prev.disabled = page === 0;
  next.disabled = page + 1 >= pages;
  prev.onclick = () => callback(page - 1);
  next.onclick = () => callback(page + 1);
  target.append(
    node(
      "span",
      `${fmt(page * size + 1)}–${fmt(Math.min(count, (page + 1) * size))} of ${fmt(count)}`,
    ),
    prev,
    next,
  );
}
function renderParsed() {
  const errors = parsed.entries.filter((e) => e.error).length,
    duplicates = parsed.entries.filter((e) => e.duplicate).length;
  $("parseSummary").textContent =
    `${fmt(parsed.variants.length)} unique · ${errors} errors · ${duplicates} repeated entries`;
  const rows = $("parseResults");
  rows.hidden = false;
  rows.replaceChildren();
  // Errors appear first so every blocked line is reachable without scrolling past good entries.
  const entries = [...parsed.entries].sort(
    (a, b) => Number(!!b.error) - Number(!!a.error) || a.line - b.line,
  );
  for (const e of entries.slice(parsePage * 50, (parsePage + 1) * 50)) {
    const row = node(
      "div",
      undefined,
      "parse-row" + (e.error ? " parse-error" : ""),
    );
    const detail = node("div");
    detail.append(
      node(
        "span",
        e.error || `${e.label}${e.duplicate ? " · repeated allele" : ""}`,
      ),
    );
    if (e.error) detail.append(node("code", e.input));
    row.append(node("span", `Line ${e.line}`), detail);
    rows.append(row);
  }
  pager("parsePager", entries.length, parsePage, 50, (p) => {
    parsePage = p;
    renderParsed();
  });
}
function checkVariants() {
  clearResults();
  invalidateParsed();
  startWorker();
  busy = true;
  controls();
  status("Parsing your variant list locally…");
  armWatchdog();
  worker.postMessage({
    type: "parse",
    text: $("variants").value,
    assembly: $("assembly").value,
  });
}
function options() {
  return {
    minVaf: Number($("minVaf").value) / 100,
    sample: $("sample").value,
    mapq: Number($("mapq").value),
    baseq: Number($("baseq").value),
    duplicates: $("duplicates").checked,
    helpers: $("helpers").checked,
    helperSupport: Number($("helperSupport").value),
    helperFraction: Number($("helperFraction").value) / 100,
    helperGap: Number($("helperGap").value),
    linkAgreement: Number($("linkAgreement").value) / 100,
  };
}
function indexComponents() {
  components = new Map();
  let component = 0;
  for (const start of graph.keys()) {
    if (components.has(start)) continue;
    const block = { id: component++, inconsistent: false, requested: 0 },
      queue = [start];
    components.set(start, { block, sign: 1 });
    for (let i = 0; i < queue.length; i++) {
      const a = queue[i],
        state = components.get(a);
      if (a < result.requestedCount) block.requested++;
      for (const edge of graph.get(a) || []) {
        const sign = state.sign * edge.sign,
          prior = components.get(edge.to);
        if (prior) {
          if (prior.sign !== sign) block.inconsistent = true;
        } else {
          components.set(edge.to, { block, sign });
          queue.push(edge.to);
        }
      }
    }
  }
}
function connection(a, b) {
  const x = components.get(a),
    y = components.get(b);
  return x && y && x.block === y.block
    ? x.block.inconsistent
      ? "Conflicting helper block"
      : "Helper link available"
    : "Unresolved";
}
function renderResult() {
  selected = null;
  for (const id of ["summary", "pairDetail", "pathDetail", "evidenceRows"])
    $(id).replaceChildren();
  $("pathDetail").hidden = $("evidence").hidden = true;
  $("results").hidden = false;
  pairPage = helperPage = 0;
  usedHelpers.clear();
  pairMap = new Map(result.pairs.map((p) => [`${p.a}:${p.b}`, p]));
  requestedPairs = result.pairs.filter(
    (p) => p.a < result.requestedCount && p.b < result.requestedCount,
  );
  totalAll = (result.requestedCount * (result.requestedCount - 1)) / 2;
  graph = helperGraph(result.variants, result.helperPairs, result.options);
  indexComponents();
  const covered = requestedPairs.filter((p) => pairMetrics(p).valid > 0).length;
  for (const [value, label, cls] of [
    [result.requestedCount, "requested variants"],
    [result.fragments, "read pairs / fragments"],
    [covered, `pairs callable at both sites / ${fmt(totalAll)} possible`],
    [result.helperCount, "qualifying helper SNPs", "helper"],
  ]) {
    const tile = node("div", undefined, "stat" + (cls ? " " + cls : ""));
    tile.append(node("strong", fmt(value)), node("span", label));
    $("summary").append(tile);
  }
  $("runSettings").textContent =
    `Sample: ${result.options.sample} · MAPQ ≥${result.options.mapq} · BQ ≥${result.options.baseq} · minimum VAF ${100 * (result.options.minVaf ?? 0.02)}% · duplicates ${result.options.duplicates ? "included in direct counts" : "excluded"} · ${fmt(result.filtered)} filtered alignments`;
  for (const id of ["variantA", "variantB"]) {
    const select = $(id);
    select.replaceChildren();
    for (const [name, helper] of [
      ["Requested variants", false],
      ["Helper SNPs", true],
    ]) {
      const group = document.createElement("optgroup");
      group.label = name;
      for (const v of result.variants.filter((v) => !!v.helper === helper)) {
        const o = node("option", `${v.helper ? "Helper · " : ""}${v.label}`);
        o.value = v.id;
        group.append(o);
      }
      if (group.children.length) select.append(group);
    }
  }
  $("variantA").value = "0";
  $("variantB").value = String(Math.min(1, result.variants.length - 1));
  $("helpersPanel").hidden = !result.helperCount;
  renderPairs();
  renderHelpers();
  if (result.requestedCount >= 2) inspect(0, 1);
  else if (!result.requestedCount)
    $("pairDetail").append(
      node(
        "p",
        "No requested variants remain above the VAF threshold. Lower the minimum VAF and rerun to include suppressed alleles.",
        "small muted",
      ),
    );
  else
    $("pairDetail").append(
      node(
        "p",
        "One requested site was assessed. Add another to compare pairwise phase; site coverage is reported below.",
        "small muted",
      ),
      node("p", siteSummary(result.sites[0]), "small"),
    );
}
function siteSummary(s) {
  return `ALT ${fmt(s.alt)} · REF ${fmt(s.ref)} · other ${fmt(s.other)} · conflict ${fmt(s.conflict)} · callable ALT fraction ${percent(s.alt + s.ref ? s.alt / (s.alt + s.ref) : null)}`;
}
function zero(a, b) {
  return {
    a,
    b,
    both: 0,
    leftOnly: 0,
    rightOnly: 0,
    neither: 0,
    other: 0,
    conflict: 0,
  };
}
function getPair(a, b) {
  return (
    pairMap.get(`${Math.min(a, b)}:${Math.max(a, b)}`) ||
    zero(Math.min(a, b), Math.max(a, b))
  );
}
function pairAt(index) {
  let a = 0,
    remaining = index,
    row = result.requestedCount - 1;
  while (remaining >= row && row > 0) {
    remaining -= row;
    a++;
    row--;
  }
  return getPair(a, a + 1 + remaining);
}
function freqCell(value, n, d) {
  const cell = node("td");
  cell.append(
    node("span", percent(value), "frequency"),
    node("small", d ? `${fmt(n)} / ${fmt(d)}` : "No ALT-bearing joint calls"),
  );
  return cell;
}
function pairName(v) {
  return node("code", v.label, v.helper ? "helper" : undefined);
}
function renderPairs() {
  const all = $("pairMode").value === "all",
    count = all ? totalAll : requestedPairs.length,
    body = $("pairRows");
  body.replaceChildren();
  if (!count) {
    const row = node("tr"),
      cell = node(
        "td",
        "No requested pair has joint observations. Switch to all pairs or inspect a pair to check helper links.",
      );
    cell.colSpan = 9;
    row.append(cell);
    body.append(row);
  }
  for (let i = pairPage * 50; i < Math.min(count, (pairPage + 1) * 50); i++) {
    const p = all ? pairAt(i) : requestedPairs[i],
      m = pairMetrics(p),
      row = document.createElement("tr");
    row.dataset.pair = `${p.a}:${p.b}`;
    const names = node("td", undefined, "pair-names");
    for (const id of [p.a, p.b]) {
      const div = node("div");
      div.append(pairName(result.variants[id]));
      names.append(div);
    }
    const evidence = node("td"),
      label = m.valid ? "Direct read pairs" : connection(p.a, p.b);
    evidence.append(
      node(
        "span",
        label,
        `tag ${label.includes("helper") || label.includes("Helper") ? "helper" : m.valid ? "" : "unknown"}`,
      ),
    );
    evidence.append(
      node("small", `${fmt(m.valid)} REF/ALT-callable fragments`),
    );
    row.append(
      names,
      freqCell(m.together, p.both, m.atLeastOne),
      freqCell(m.separate, p.leftOnly + p.rightOnly, m.atLeastOne),
    );
    for (const n of [
      p.both,
      p.leftOnly,
      p.rightOnly,
      p.neither,
      p.other + p.conflict,
    ])
      row.append(node("td", fmt(n)));
    row.append(evidence);
    row.onclick = () => inspect(p.a, p.b);
    body.append(row);
  }
  pager("pairPager", count, pairPage, 50, (p) => {
    pairPage = p;
    renderPairs();
  });
}
function renderHelpers() {
  const helpers = result.variants.filter((v) => v.helper),
    body = $("helperRows");
  body.replaceChildren();
  for (const v of helpers.slice(helperPage * 50, (helperPage + 1) * 50)) {
    const row = node(
        "tr",
        undefined,
        usedHelpers.has(v.id) ? "used-helper" : undefined,
      ),
      label = node("td");
    label.append(pairName(v));
    row.append(
      label,
      node("td", fmt(v.altCount)),
      node("td", fmt(v.depth)),
      node("td", percent(v.fraction)),
      node(
        "td",
        usedHelpers.has(v.id) ? "Used in this path" : "Not in selected path",
      ),
    );
    body.append(row);
  }
  pager("helperPager", helpers.length, helperPage, 50, (p) => {
    helperPage = p;
    renderHelpers();
  });
}
function inspect(a, b, keepPath = false) {
  a = Number(a);
  b = Number(b);
  if (a === b) {
    $("pairDetail").replaceChildren(
      node("p", "Choose two different variants.", "small"),
    );
    $("pathDetail").hidden = $("evidence").hidden = true;
    selected = null;
    return;
  }
  if (a > b) [a, b] = [b, a];
  selected = [a, b];
  $("variantA").value = a;
  $("variantB").value = b;
  const p = getPair(a, b),
    m = pairMetrics(p),
    detail = $("pairDetail");
  detail.replaceChildren();
  const names = node("p");
  names.append(
    pairName(result.variants[a]),
    document.createTextNode(" ↔ "),
    pairName(result.variants[b]),
  );
  detail.append(names);
  const grid = node("div", undefined, "detail-grid");
  for (const [label, value, n, d] of [
    ["Together", m.together, p.both, m.atLeastOne],
    ["Separate", m.separate, p.leftOnly + p.rightOnly, m.atLeastOne],
    ["B ALT among A ALT", m.bGivenA, p.both, p.both + p.leftOnly],
    ["A ALT among B ALT", m.aGivenB, p.both, p.both + p.rightOnly],
  ]) {
    const tile = node("div");
    tile.append(
      node("strong", percent(value)),
      node("span", label),
      node("small", `${n} / ${d} fragments`),
    );
    grid.append(tile);
  }
  detail.append(
    grid,
    node(
      "p",
      `Direct calls: both ALT ${p.both}, A-only ${p.leftOnly}, B-only ${p.rightOnly}, neither ${p.neither}. Excluded: other allele ${p.other}, conflicting mates ${p.conflict}. ${m.valid ? "" : "No directly measured joint REF/ALT observations."}`,
      "small",
    ),
  );
  for (const id of [a, b])
    detail.append(
      node(
        "p",
        `${id === a ? "A" : "B"} site · ${siteSummary(result.sites[id])}`,
        "small muted",
      ),
    );
  if (result.variants[a].chr !== result.variants[b].chr)
    detail.append(
      node(
        "p",
        "Different contigs: no local cis frequency is reported, even if a discordant pair maps to both.",
        "small",
      ),
    );
  if (!keepPath) renderPath(findHelperPath(graph, result.variants, a, b));
  $("evidenceRows").replaceChildren();
  $("evidencePager").replaceChildren();
  $("evidence").hidden = false;
  worker.postMessage({ type: "evidence", a, b, offset: 0 });
}
function renderPath(path) {
  const target = $("pathDetail");
  target.replaceChildren();
  target.hidden = !result.options.helpers;
  usedHelpers = new Set(path.helperIds || []);
  renderHelpers();
  if (!result.options.helpers) return;
  target.append(
    node(
      "h4",
      path.status === "linked"
        ? "Helper-assisted phase link"
        : path.status === "conflicting"
          ? "Conflicting helper block"
          : "No qualifying helper path",
    ),
  );
  if (path.status === "unresolved") {
    target.append(
      node(
        "p",
        "No additional indirect link was established through qualifying helper SNPs. Nearby markers, coverage gaps or insufficient link support can leave this pair unresolved.",
        "small",
      ),
    );
    return;
  }
  if (path.status === "conflicting")
    target.append(
      node(
        "p",
        "The reachable marker block contains contradictory phase paths. This chain is shown for inspection; it is not a resolved phase assignment.",
        "small",
      ),
    );
  else
    target.append(
      node(
        "p",
        `${path.sign === 1 ? "Coupled ALT-marker relationship" : "Opposite ALT-marker relationship"} along ${path.edges.length} local links. This is indirect evidence; an end-to-end same-allele frequency cannot be measured from this chain.`,
        "small",
      ),
    );
  const nodes = node("div", undefined, "path-nodes");
  for (const [i, id] of path.nodes.entries()) {
    if (i) nodes.append(node("span", "→"));
    const v = result.variants[id];
    nodes.append(
      node(
        "span",
        `${v.helper ? "Helper · " : ""}${v.label}`,
        "path-node" + (v.helper ? " helper" : ""),
      ),
    );
  }
  target.append(nodes);
  const edges = node("div", undefined, "path-edges");
  for (const [i, edge] of path.edges.entries()) {
    const button = node(
      "button",
      `Link ${i + 1}: ${edge.sign === 1 ? "coupled" : "opposite"} · ${edge.count}/${edge.total} (${percent(edge.agreement)})`,
    );
    button.onclick = () => inspect(edge.a, edge.b, true);
    edges.append(button);
  }
  const back = node("button", "Return to requested pair");
  back.onclick = () => inspect(path.nodes[0], path.nodes.at(-1));
  edges.append(back);
  target.append(
    edges,
    node(
      "p",
      `Purple markers were used in this path. Each link uses coordinate/UMI-distinct molecules, ≥${result.options.helperSupport} observations of each supported phase pattern, and ≥${percent(result.options.linkAgreement)} agreement. Counts below show raw read-pair templates; they can exceed the deduplicated link counts.`,
      "small",
    ),
  );
}
function pairedTrack(fragment, a, b, domain) {
  if (a.chr !== b.chr || fragment.ends.some((e) => e.chr !== a.chr))
    return null;
  const canvas = document.createElement("canvas");
  canvas.width = 500;
  canvas.height = 42;
  canvas.className = "pair-track";
  canvas.setAttribute("role", "img");
  canvas.setAttribute(
    "aria-label",
    "Paired alignments on a shared genomic axis; dashed insert gap, teal requested sites, purple helper sites",
  );
  const ctx = canvas.getContext("2d"),
    scale = (p) => 10 + ((p - domain[0]) * 480) / (domain[1] - domain[0]);
  const ends = [...fragment.ends].sort((x, y) => x.start - y.start);
  if (ends.length === 2) {
    ctx.strokeStyle = "#a5b8bd";
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(scale(ends[0].end), 27);
    ctx.lineTo(scale(ends[1].start), 27);
    ctx.stroke();
    ctx.setLineDash([]);
  }
  for (const e of ends) {
    ctx.fillStyle = e.endFlag === 128 ? "#83aeb3" : "#4c808b";
    ctx.fillRect(
      scale(e.start),
      20,
      Math.max(2, scale(e.end) - scale(e.start)),
      13,
    );
  }
  for (const [v, label] of [
    [a, "A"],
    [b, "B"],
  ]) {
    ctx.fillStyle = ctx.strokeStyle = v.helper ? "#7944ab" : "#116c65";
    const x = scale(v.pos + 0.5);
    ctx.beginPath();
    ctx.moveTo(x, 15);
    ctx.lineTo(x, 38);
    ctx.stroke();
    ctx.font = "bold 10px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(label, x, 11);
  }
  return canvas;
}
function renderEvidence(m) {
  const a = result.variants[m.a],
    b = result.variants[m.b];
  const points = [
    a.pos,
    b.pos,
    ...m.rows.flatMap((f) => f.ends.flatMap((e) => [e.start, e.end])),
  ];
  const domain = [Math.min(...points) - 5, Math.max(...points) + 5];
  const body = $("evidenceRows");
  body.replaceChildren();
  if (!m.rows.length) {
    const row = node("tr"),
      cell = node(
        "td",
        "No fragments with calls at both sites. Uncovered sites and quality-filtered calls do not become reference observations.",
      );
    cell.colSpan = 4;
    row.append(cell);
    body.append(row);
  }
  for (const f of m.rows) {
    const row = node("tr"),
      name = node("td", f.name, "read-name"),
      ends = node("td");
    name.append(node("small", `RG: ${f.rg || "none"} · ${f.sample}`));
    const track = pairedTrack(f, a, b, domain);
    if (track) ends.append(track);
    for (const e of [...f.ends].sort((a, b) => a.endFlag - b.endFlag)) {
      const label = e.endFlag === 64 ? "R1" : e.endFlag === 128 ? "R2" : "Read";
      ends.append(
        node(
          "div",
          `${label} ${e.chr}:${e.start + 1}–${e.end} (${e.strand}) · MQ ${e.mq} · ${e.cigar}`,
          "align-end",
        ),
      );
    }
    if (f.ends.length === 1)
      ends.append(
        node("small", "Only one end contributes within the queried scope"),
      );
    row.append(name, ends);
    for (const call of f.calls) {
      const td = node("td");
      td.append(
        node(
          "span",
          { A: "ALT", R: "REF", O: "OTHER", X: "CONFLICT" }[call] || "UNKNOWN",
          `call ${call}`,
        ),
      );
      row.append(td);
    }
    body.append(row);
  }
  pager("evidencePager", m.total, Math.floor(m.offset / 50), 50, (p) =>
    worker.postMessage({ type: "evidence", a: m.a, b: m.b, offset: p * 50 }),
  );
}
$("open").onclick = () => {
  stopWorker();
  opened = false;
  clearResults();
  invalidateParsed();
  startWorker();
  busy = true;
  controls();
  status("Opening local BAM header and index…");
  armWatchdog();
  worker.postMessage({
    type: "open",
    files: { bam: $("bam").files[0], index: $("index").files[0] },
  });
};
$("parse").onclick = checkVariants;
$("run").onclick = () => {
  for (const id of [
    "minVaf",
    "mapq",
    "baseq",
    "helperSupport",
    "helperFraction",
    "helperGap",
    "linkAgreement",
  ])
    if (!$(id).checkValidity()) {
      $(id).reportValidity();
      return;
    }
  clearResults();
  busy = true;
  controls();
  status("Reading requested loci locally…");
  armWatchdog();
  worker.postMessage({ type: "run", options: options() });
};
$("clear").onclick = reset;
$("cancel").onclick = () => {
  stopWorker();
  opened = false;
  invalidateParsed();
  clearResults();
  finish();
  $("fileStatus").textContent = "Processing canceled; reopen files";
  status(
    "Processing canceled and analysis memory released. Reopen files to start again.",
  );
};
$("inspect").onclick = () => inspect($("variantA").value, $("variantB").value);
$("pairMode").onchange = () => {
  pairPage = 0;
  if (result) renderPairs();
};
for (const id of ["variants", "assembly"])
  $(id).addEventListener(id === "variants" ? "input" : "change", () => {
    clearResults();
    invalidateParsed();
    controls();
    status("Variant list changed. Check it before running.");
  });
for (const id of ["bam", "index"])
  $(id).addEventListener("change", () => {
    stopWorker();
    opened = false;
    clearResults();
    invalidateParsed();
    $("fileStatus").textContent = "Files changed; open them to continue";
    controls();
  });
for (const id of [
  "sample",
  "mapq",
  "baseq",
  "duplicates",
  "helpers",
  "helperSupport",
  "helperFraction",
  "helperGap",
  "linkAgreement",
])
  $(id).addEventListener("change", () => {
    clearResults();
    controls();
    status("Evidence settings changed. Run phasing to update results.");
  });
watchVaf($("minVaf"), {
  current: () => result?.options.minVaf,
  busy: () => busy,
  update: (minVaf) => {
    busy = true;
    controls();
    armWatchdog();
    status("Updating visible variants from loaded reads…");
    worker.postMessage({ type: "vaf", minVaf });
  },
  error: (e) => status(e.message, true),
});
window.addEventListener("pagehide", reset);
reset();
status("Open BAM/BAI and check your variant list to begin.");
