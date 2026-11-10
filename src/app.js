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
import { cigarOps } from "./analysis.js";
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
      data = m;
      view = { start: m.region.start, end: m.region.end };
      rows = 0;
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
  worker?.terminate();
  worker = null;
  if (url) URL.revokeObjectURL(url);
  url = null;
  clearTimeout(timer);
  opened = false;
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
  for (const id of ["locus", "transcript", "sample", "v1", "v2"])
    $(id).replaceChildren();
  for (const id of [
    "blocks",
    "variants",
    "events",
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
  const r = data.result,
    group = new Map();
  for (const f of r.fragments)
    for (const id of f.reads)
      group.set(id, f.groups.join(", ") || "Unassigned");
  let reads = r.reads
    .filter((x) => x.start < view.end && x.end > view.start)
    .filter(
      (x) =>
        $("layer").value === "all" ||
        ($("layer").value === "excluded" ? x.excluded : !x.excluded),
    );
  const sort = $("sort").value;
  reads.sort((a, b) =>
    sort === "group"
      ? (group.get(a.id) || "Excluded").localeCompare(
          group.get(b.id) || "Excluded",
        ) || a.start - b.start
      : sort === "quality"
        ? b.mq - a.mq
        : sort === "strand"
          ? (a.flags & 16) - (b.flags & 16) || a.start - b.start
          : a.start - b.start,
  );
  return reads.map((x) => ({ ...x, group: group.get(x.id) || "Excluded" }));
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
        ? "Overlapping allele · inspect locally"
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
  draw();
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
  const canvas = $("pileup"),
    width = Math.max(400, $("canvasWrap").clientWidth),
    all = visibleReads();
  rows = Math.min(rows, Math.max(0, Math.ceil(all.length / 150) - 1) * 150);
  painted = all.slice(rows, rows + 150);
  const height = 110 + Math.max(5, painted.length) * 15,
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
  const max = Math.max(1, ...data.result.coverage);
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
  painted.forEach((r, i) => {
    const y = 110 + i * 15;
    c.fillStyle = r.excluded ? "#a68d80" : "#638975";
    c.fillText(r.group.slice(0, 16), 8, y + 9);
    let p = r.start,
      q = 0;
    for (const [n, op] of cigarOps(r.cigar)) {
      if ("M=X".includes(op)) {
        c.fillStyle = r.excluded ? "#dfe2de" : "#c6d0c8";
        c.fillRect(
          Math.max(left, x(p)),
          y,
          Math.max(0, Math.min(width - 12, x(p + n)) - Math.max(left, x(p))),
          10,
        );
        for (let j = 0; j < n; j++) {
          const pos = p + j;
          if (pos < view.start || pos >= view.end) continue;
          const b = r.seq[q + j];
          if (b !== data.reference[pos - data.region.start]) {
            c.fillStyle = colors[b] || "#aaa";
            c.fillRect(x(pos), y, Math.max(1, scale), 10);
          }
          if (scale > 9) {
            c.fillStyle =
              b === data.reference[pos - data.region.start]
                ? "#5f7065"
                : "#fff";
            c.fillText(b, x(pos), y + 9);
          }
        }
        p += n;
        q += n;
      } else if (op === "D" || op === "N") {
        c.strokeStyle = op === "D" ? "#7b817c" : "#bcc9c0";
        c.beginPath();
        c.moveTo(x(p), y + 5);
        c.lineTo(x(p + n), y + 5);
        c.stroke();
        p += n;
      } else if (op === "I") {
        c.fillStyle = "#9c70b3";
        c.fillRect(x(p) - 2, y - 2, 4, 14);
        q += n;
      } else if (op === "S") {
        c.fillStyle = "#b3a0c2";
        const xp = x(p);
        c.fillRect(xp - 2, y, 4, 10);
        if ($("clips").checked && scale > 8) {
          c.fillText(r.seq.slice(q, q + n), xp, y + 10);
        }
        q += n;
      }
    }
  });
  $("rowCount").textContent =
    `Records ${all.length ? rows + 1 : 0}–${Math.min(rows + 150, all.length)} of ${all.length} in view · ${Math.round(view.start + 1)}–${Math.round(view.end)}`;
  $("readTable").replaceChildren();
  for (const r of painted) {
    const tr = el("tr");
    for (const v of [
      `Record ${r.id + 1}`,
      r.start + 1,
      r.cigar,
      r.mq,
      r.flags,
      r.group,
    ])
      tr.append(el("td", v));
    tr.tabIndex = 0;
    const detail = () => {
      $("readDetail").textContent =
        `Read name (selected explicitly): ${r.name}\nStart: ${r.start + 1}\nCIGAR: ${r.cigar}\nFlags: ${r.flags}\nMate start: ${r.mate >= 0 ? r.mate + 1 : "unavailable"}\nSequence: ${r.seq}\nBase qualities: ${r.qual.join(",")}`;
    };
    tr.onclick = detail;
    tr.onkeydown = (e) => {
      if (e.key === "Enter") detail();
    };
    $("readTable").append(tr);
  }
}
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
    draw();
  }
};
for (const id of ["sort", "layer", "clips"])
  $(id).onchange = () => {
    rows = 0;
    draw();
  };
$("prevRows").onclick = () => {
  rows = Math.max(0, rows - 150);
  draw();
};
$("nextRows").onclick = () => {
  rows += 150;
  draw();
};
$("pileup").onkeydown = (e) => {
  if (["ArrowLeft", "ArrowRight", "+", "-"].includes(e.key)) {
    e.preventDefault();
    if (e.key === "ArrowLeft") pan(-0.25);
    else if (e.key === "ArrowRight") pan(0.25);
    else zoom(e.key === "+" ? 0.5 : 2);
  }
};
new ResizeObserver(() => draw()).observe($("canvasWrap"));
