import { test } from "node:test";
import assert from "node:assert/strict";
import {
  IGV_COLORS,
  IGV_INSERTION,
  IGV_SKIP,
  qualityAlpha,
  drawIGVBase,
  drawIGVInsertion,
  drawIGVGap,
} from "../v3/igv-render.js";
const ctx = () => ({
  globalAlpha: 1,
  ops: [],
  fillRect(...v) {
    this.ops.push(["rect", this.fillStyle, this.globalAlpha, ...v]);
  },
  fillText(...v) {
    this.ops.push(["text", this.fillStyle, this.globalAlpha, ...v]);
  },
  measureText(t) {
    return { width: t.length * 6 };
  },
  setLineDash(v) {
    this.ops.push(["dash", v]);
  },
  beginPath() {},
  moveTo() {},
  lineTo() {},
  stroke() {
    this.ops.push(["stroke", this.strokeStyle, this.lineWidth]);
  },
});
test("desktop IGV shipped light-theme palette is exact, including orange G", () =>
  assert.deepEqual(IGV_COLORS, {
    A: "#009600",
    C: "#0000ff",
    G: "#d17105",
    T: "#ff0000",
    N: "#808080",
  }));
test("soft-clip bases use each nucleotide color as blocks at overview zoom and colored letters at base zoom", () => {
  const c = ctx();
  for (const [i, b] of [..."TGCA"].entries())
    drawIGVBase(c, b, i * 3, 0, 3, 18, 35);
  assert.deepEqual(
    c.ops.map((o) => o[1]),
    ["#ff0000", "#d17105", "#0000ff", "#009600"],
  );
  const z = ctx();
  drawIGVBase(z, "A", 0, 0, 12, 18, 35);
  assert.equal(z.ops[0][0], "text");
  assert.equal(z.ops[0][1], "#009600");
});
test("IGV base-quality shading is quantized to tenths with Q5/Q20 defaults", () => {
  assert.equal(qualityAlpha(0), 0.2);
  assert.equal(qualityAlpha(5), 0.2);
  assert.equal(qualityAlpha(10), 0.4);
  assert.equal(qualityAlpha(20), 1);
  const c = ctx();
  drawIGVBase(c, "T", 0, 0, 2, 18, 10);
  assert.equal(c.ops[0][2], 0.4);
  assert.equal(c.globalAlpha, 1);
});
test("insertion marker is a purple I with wings, not a triangle; size labels appear only when they fit", () => {
  const c = ctx();
  drawIGVInsertion(c, 50, 0, 18, 1, 2);
  assert.equal(c.ops.filter((o) => o[0] === "rect").length, 3);
  assert.ok(c.ops.every((o) => o[1] === IGV_INSERTION));
  const z = ctx();
  drawIGVInsertion(z, 50, 0, 18, 10, 8);
  assert.ok(
    z.ops.some((o) => o[0] === "text" && o[1] === "#ffffff" && o[3] === "10"),
  );
});
test("deletions are black with a centered white size label and skipped regions use the desktop blue-gray", () => {
  const c = ctx();
  drawIGVGap(c, 0, 100, 0, 18, 8, "D");
  assert.deepEqual(
    c.ops.find((o) => o[0] === "stroke"),
    ["stroke", "#000000", 2],
  );
  assert.ok(
    c.ops.some(
      (o) => o[0] === "text" && o[1] === IGV_INSERTION && o[3] === "8",
    ),
  );
  const n = ctx();
  drawIGVGap(n, 0, 100, 0, 18, 100, "N");
  assert.equal(n.ops.find((o) => o[0] === "stroke")[1], IGV_SKIP);
  assert.ok(!n.ops.some((o) => o[0] === "text"));
});
