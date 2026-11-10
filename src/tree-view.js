import { cutHierarchy } from "./hierarchy.js";
const palette = [
  "#277d75",
  "#8861a2",
  "#bd7846",
  "#497ab0",
  "#a65869",
  "#718c42",
  "#82714e",
  "#547e91",
];
export const groupColor = (block, group) =>
  palette[(block * 3 + group - 1) % palette.length];
const pct = (n) => (n === null ? "unknown" : (100 * n).toFixed(1) + "%");
export function setupHierarchy({ getData, onSelect }) {
  const $ = (id) => document.getElementById(id);
  let hits = [],
    selected = null,
    lastHover = "";
  const block = () =>
    getData()?.result.hierarchy?.find(
      (b) => b.blockID === +$("treeBlock").value,
    );
  function info(b, g) {
    const comparisons = b.groups
      .filter((h) => h.id !== g.id)
      .map((h) => {
        const p = b.pairs[b.groups.indexOf(g)][b.groups.indexOf(h)];
        return p.distance === null
          ? `H${h.id}: no jointly observed SNPs`
          : `H${h.id}: ${p.sharedAlternate} shared ALT, ${p.sharedReference} shared REF, ${p.different}/${p.compared} SNPs differ (${pct(p.distance)})`;
      });
    return `B${b.blockID} H${g.id}\n${g.fragments} assigned eligible fragments / ${b.informative} informative (${pct(g.assignedShare)}); ${g.records} alignment records.\n${g.ambiguousCompatible} additional fragments compatible but ambiguous (not assigned).\nObserved site allele fractions across all eligible informative fragments:\n${g.alleleFractions.length ? g.alleleFractions.map((v) => `${v.pos + 1} ${v.ref || "∅"}→${v.allele}: ${v.alt}/${v.total} (${pct(v.fraction)})`).join("\n") : "No alternate allele in this pattern."}\n${comparisons.length ? "SNP similarity to other local groups:\n" + comparisons.join("\n") : "Only one supported local group."}\nGroup share is not a cancer-cell fraction; site allele fractions can be shared by several groups.`;
  }
  function showGroup(blockID, groupID) {
    const b = getData()?.result.hierarchy?.find((b) => b.blockID === blockID),
      g = b?.groups.find((g) => g.id === groupID);
    if (!g) return;
    const key = blockID + ":" + groupID;
    if (lastHover === key) return;
    lastHover = key;
    selected = { blockID, groupID };
    $("treeBlock").value = String(blockID);
    $("groupInfo").textContent = info(b, g);
    draw();
  }
  function draw() {
    const b = block(),
      canvas = $("treeCanvas"),
      width = Math.max(
        210,
        Math.min(600, $("hierarchyPanel").clientWidth - 30),
      ),
      height = b ? Math.max(80, b.groups.length * 44 + 38) : 80,
      dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    canvas.style.height = height + "px";
    const c = canvas.getContext("2d");
    c.scale(dpr, dpr);
    c.clearRect(0, 0, width, height);
    hits = [];
    const cut = +$("treeCut").value / 100;
    $("cutLabel").textContent = `${Math.round(cut * 100)}% SNP difference`;
    if (!b) {
      $("treeSummary").textContent =
        "No phase groups to cluster in this interval.";
      return;
    }
    const clusters = cutHierarchy(b.roots, cut);
    $("treeSummary").textContent =
      `B${b.blockID}: ${b.groups.length} supported groups → ${clusters.length} similarity clusters at this cut. ${b.sites.length} distinct SNP sites. ${b.ambiguous} ambiguous; ${b.conflicting} unresolved fragments. Fractions use ${b.informative} informative fragments.`;
    $("blockChips").replaceChildren();
    for (const entry of getData().result.hierarchy) {
      const button = document.createElement("button");
      button.textContent = `B${entry.blockID}: ${entry.groups.length} groups`;
      button.className = entry.blockID === b.blockID ? "selected" : "";
      button.onclick = () => {
        $("treeBlock").value = String(entry.blockID);
        lastHover = "";
        selected = null;
        $("groupInfo").textContent =
          "Hover a group for counts and SNP comparisons.";
        draw();
      };
      $("blockChips").append(button);
    }
    const treeWidth = Math.max(55, width * 0.36),
      left = 8,
      leafX = left + treeWidth,
      x = (h) => leafX - treeWidth * h;
    let leaf = 0;
    const positions = new Map();
    function layout(n) {
      if (!n.left) {
        const y = 35 + leaf++ * 44;
        positions.set(n.id, y);
        return y;
      }
      const y = (layout(n.left) + layout(n.right)) / 2;
      positions.set(n.id, y);
      return y;
    }
    b.roots.forEach(layout);
    c.font = "9px sans-serif";
    c.fillStyle = "#789081";
    c.fillText("100% different", left, 11);
    c.fillText("0%", leafX - 8, 11);
    c.strokeStyle = "#b1c4b7";
    c.setLineDash([3, 3]);
    c.beginPath();
    c.moveTo(x(cut), 17);
    c.lineTo(x(cut), height);
    c.stroke();
    c.setLineDash([]);
    function paint(n) {
      const y = positions.get(n.id);
      if (n.left) {
        const a = positions.get(n.left.id),
          d = positions.get(n.right.id);
        c.strokeStyle = n.height <= cut ? "#4f8971" : "#9ca9a0";
        c.beginPath();
        c.moveTo(x(n.left.height), a);
        c.lineTo(x(n.height), a);
        c.lineTo(x(n.height), d);
        c.lineTo(x(n.right.height), d);
        c.stroke();
        paint(n.left);
        paint(n.right);
      } else {
        const g = b.groups.find((g) => g.id === n.groupID),
          active = selected?.blockID === b.blockID && selected.groupID === g.id,
          color = groupColor(b.blockID, g.id);
        if (active) {
          c.fillStyle = "#e7f2e9";
          c.fillRect(0, y - 13, width, 40);
        }
        c.fillStyle = color;
        c.beginPath();
        c.arc(leafX, y, 4, 0, Math.PI * 2);
        c.fill();
        c.font = "11px sans-serif";
        c.fillText(`H${g.id} · ${g.fragments} fragments`, leafX + 9, y + 3);
        const fractions = g.alleleFractions
          .map((v) => v.fraction)
          .filter((v) => v !== null);
        const range = fractions.length
          ? `${Math.round(Math.min(...fractions) * 100)}–${Math.round(Math.max(...fractions) * 100)}%`
          : "none";
        c.font = "9px sans-serif";
        c.fillStyle = "#657e6e";
        c.fillText(
          `${pct(g.assignedShare)} group share · site VAF ${range}`,
          leafX + 9,
          y + 16,
          width - leafX - 11,
        );
        c.fillStyle = color;
        c.fillRect(
          leafX + 9,
          y + 21,
          (width - leafX - 15) * (g.assignedShare || 0),
          3,
        );
        hits.push({ y, group: g });
      }
    }
    b.roots.forEach(paint);
  }
  function render() {
    lastHover = "";
    selected = null;
    $("groupInfo").textContent =
      "Hover a group in either graph for counts, allele fractions, and SNP comparisons.";
    $("treeBlock").replaceChildren();
    for (const b of getData()?.result.hierarchy || [])
      $("treeBlock").append(
        new Option(
          `B${b.blockID} · ${b.groups.length} groups · ${b.start + 1}–${b.end}`,
          b.blockID,
        ),
      );
    draw();
  }
  $("treeBlock").onchange = () => {
    lastHover = "";
    selected = null;
    $("groupInfo").textContent =
      "Hover a group for counts and SNP comparisons.";
    draw();
  };
  $("treeCut").oninput = draw;
  $("treeCanvas").onpointermove = (e) => {
    const y = e.clientY - $("treeCanvas").getBoundingClientRect().top,
      hit = hits.find((h) => Math.abs(h.y - y) < 22),
      b = block();
    if (hit && b) {
      showGroup(b.blockID, hit.group.id);
      $("treeCanvas").title = info(b, hit.group);
    }
  };
  $("treeCanvas").onclick = (e) => {
    const y = e.clientY - $("treeCanvas").getBoundingClientRect().top,
      hit = hits.find((h) => Math.abs(h.y - y) < 22),
      b = block();
    if (hit && b) {
      showGroup(b.blockID, hit.group.id);
      onSelect(b.blockID, hit.group.id);
    }
  };
  $("treeCanvas").onkeydown = (e) => {
    const b = block();
    if (!b?.groups.length) return;
    if (["ArrowDown", "ArrowUp", "Enter"].includes(e.key)) {
      e.preventDefault();
      const i = b.groups.findIndex((g) => g.id === selected?.groupID),
        next =
          e.key === "ArrowUp"
            ? (i - 1 + b.groups.length) % b.groups.length
            : e.key === "ArrowDown"
              ? (i + 1) % b.groups.length
              : Math.max(i, 0),
        g = b.groups[next];
      showGroup(b.blockID, g.id);
      if (e.key === "Enter") onSelect(b.blockID, g.id);
    }
  };
  function clear() {
    hits = [];
    selected = null;
    lastHover = "";
    for (const id of ["treeBlock", "blockChips", "treeSummary", "groupInfo"])
      $(id).replaceChildren();
    $("treeCanvas").width = 1;
    $("treeCanvas").height = 1;
    $("treeCanvas").removeAttribute("title");
  }
  function showSummary(text) {
    if (lastHover === text) return;
    lastHover = text;
    selected = null;
    $("groupInfo").textContent = text;
    draw();
  }
  return { render, draw, clear, showGroup, showSummary, info, block };
}
