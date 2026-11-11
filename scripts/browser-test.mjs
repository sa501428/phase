import { chromium } from "playwright-core";
import { resolve } from "node:path";
import { writeFile, mkdir, readFile } from "node:fs/promises";
import assert from "node:assert/strict";
import { fixture } from "./fixtures.mjs";
// Deliberately separated pairs test shared-lane selection, not just deep
// coverage where every template happens to overlap every other template.
const reference = "A".repeat(800),
  reads = [];
for (const [name, start, duplicate] of [
  ["left", 1020, false],
  ["right", 1240, false],
  ["overlap", 1060, false],
  ["clipped", 1440, false],
  ["duplicate", 1020, true],
]) {
  for (const [offset, flags] of [
    [0, 99],
    [40, 147],
  ]) {
    let seq = "A".repeat(60),
      cigar = "60M";
    if (name === "left" && offset === 0)
      seq = seq.slice(0, 30) + "C" + seq.slice(31);
    if (name === "clipped" && offset === 0) {
      seq = "TGCA" + seq;
      cigar = "4S60M";
    }
    if (name === "left") {
      const at = 1070 - start - offset;
      seq = seq.slice(0, at) + "ACGT" + seq.slice(at);
      cigar = `${at}M4I${60 - at}M`;
    }
    reads.push({
      id: reads.length,
      name,
      start: start + offset,
      end: start + offset + 60,
      seq,
      cigar,
      qual: Array(seq.length).fill(35),
      flags: flags | (duplicate ? 1024 : 0),
      mq: 60,
      rg: "synthetic",
      sample: "synthetic",
      mate: start + (offset ? 0 : 40),
      mateRef: 0,
    });
  }
}
const data = fixture({
  reads,
  reference,
  region: { chr: "synthetic", start: 1000, end: 1800 },
});
await mkdir("/tmp/phase-packed-tests", { recursive: true });
for (const [ext, bytes] of Object.entries(data))
  await writeFile(
    "/tmp/phase-packed-tests/packed." +
      ({ twobit: "2bit", fasta: "fa", fai: "fa.fai" }[ext] || ext),
    bytes,
  );
const browser = await chromium.launch({
  executablePath:
    process.env.CHROME_PATH ||
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  headless: true,
  args: ["--disable-background-networking"],
});
try {
  const page = await browser.newPage({
      viewport: { width: 1600, height: 1100 },
    }),
    errors = [],
    requests = [],
    consoleErrors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  page.on("request", (r) => {
    if (/^https?:/.test(r.url())) requests.push(r.url());
  });
  page.on("console", (m) => {
    if (m.type() === "error") consoleErrors.push(m.text());
  });
  await page.addInitScript(() => {
    const textOriginal = CanvasRenderingContext2D.prototype.fillText;
    CanvasRenderingContext2D.prototype.fillText = function (
      text,
      x,
      y,
      ...args
    ) {
      if (["pileup", "referenceTrack"].includes(this.canvas.id))
        (window.__texts ??= []).push({
          id: this.canvas.id,
          text,
          x,
          y,
          align: this.textAlign,
        });
      return textOriginal.call(this, text, x, y, ...args);
    };
    const original = CanvasRenderingContext2D.prototype.fillRect;
    CanvasRenderingContext2D.prototype.fillRect = function (x, y, w, h) {
      if (this.canvas.id === "pileup") {
        if (x === 0 && y === 0 && this.fillStyle === "#ffffff")
          window.__rects = [];
        (window.__rects ??= []).push({ x, y, w, h, color: this.fillStyle });
      }
      return original.call(this, x, y, w, h);
    };
  });
  await page.goto("file://" + resolve("index.html"));
  await page.setInputFiles("#bam", "/tmp/phase-packed-tests/packed.bam");
  await page.setInputFiles("#index", "/tmp/phase-packed-tests/packed.bai");
  await page.setInputFiles("#reference", "/tmp/phase-packed-tests/packed.2bit");
  await page.click("#open");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("Opened indexed BAM"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#load");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("Loaded 5 read pairs"),
  );
  const canvas = page.locator("#pileup");
  await page.waitForFunction(
    () => document.querySelector("#pileup").dataset.laneCount === "2",
  );
  assert.match(
    await page.textContent("#rowCount"),
    /2 packed lanes · 4 visible pairs · 1 duplicate-flagged/,
  );
  assert.equal(await page.inputValue("#readSize"), "18");
  assert.equal(await page.locator("body").getAttribute("class"), "top-compact");
  const beforeWidth = (await canvas.boundingBox()).width;
  assert.ok(beforeWidth > 1400);
  await page.evaluate(() => new Promise(requestAnimationFrame));
  await page.screenshot({
    path: "/tmp/phase-packed-initial.png",
    fullPage: true,
  });
  const geometry = await page.evaluate(() =>
    Object.fromEntries(
      ["pileup", "readScroll"].map((id) => [
        id,
        document.getElementById(id).getBoundingClientRect().height,
      ]),
    ),
  );

  assert.ok(geometry.pileup > 450);
  assert.ok(Math.abs(geometry.pileup - geometry.readScroll) < 2);
  assert.equal(await page.locator("#inspector").isVisible(), false);
  assert.equal(await page.locator("#bam").isVisible(), false);
  const info = async (name) => {
    await page.locator("#accessible").evaluate((e) => (e.open = true));
    await page.evaluate(() => new Promise(requestAnimationFrame));
    const button = page
      .locator("#accessibleRows button")
      .filter({ hasText: new RegExp("^" + name + "(?: ✓)?$") });
    const result = await button.evaluate((e) => ({
      x: +e.dataset.x,
      y: +e.dataset.y,
      key: e.dataset.key,
    }));
    await page.locator("#accessible").evaluate((e) => (e.open = false));
    await page.waitForFunction(
      () => document.querySelector("#pileup").clientHeight > 100,
    );
    return result;
  };
  const point = async (name) => {
    const p = await info(name);
    if (name === "left")
      p.x = await page.evaluate(() => {
        const c = document.querySelector("#pileup"),
          left = +c.dataset.leftMargin;
        return (
          left +
          ((1040 - +c.dataset.viewStart) /
            (+c.dataset.viewEnd - +c.dataset.viewStart)) *
            (c.clientWidth - left - 12)
        );
      });
    await page.locator("#readScroll").evaluate((e, y) => {
      if (y < e.scrollTop || y > e.scrollTop + e.clientHeight - 10)
        e.scrollTop = Math.max(0, y - 30);
    }, p.y);
    const b = await canvas.boundingBox(),
      scroll = await page.locator("#readScroll").evaluate((e) => e.scrollTop);
    return { ...p, screenX: b.x + p.x, screenY: b.y + p.y - scroll };
  };
  const clipColors = await page.evaluate(() => {
    const c = document.querySelector("#pileup"),
      left = +c.dataset.leftMargin,
      span = +c.dataset.viewEnd - +c.dataset.viewStart,
      scale = (c.clientWidth - left - 12) / span;
    return [..."TGCA"].map((base, i) => ({
      base,
      expected: { T: "#ff0000", G: "#d17105", C: "#0000ff", A: "#009600" }[
        base
      ],
      x: left + (1436 + i - +c.dataset.viewStart) * scale,
      found: window.__rects
        .filter(
          (r) =>
            Math.abs(r.x - (left + (1436 + i - +c.dataset.viewStart) * scale)) <
            0.1,
        )
        .map((r) => r.color),
    }));
  });
  for (const c of clipColors)
    assert.ok(c.found.includes(c.expected), JSON.stringify(c));
  const left = await point("left"),
    right = await point("right"),
    overlap = await point("overlap");
  assert.equal(left.y, right.y);
  assert.notEqual(left.y, overlap.y);
  // Plain movement never opens details by default. Right-click first offers a
  // menu; details and optional hover behavior are explicit user actions.
  await page.mouse.move(left.screenX, left.screenY);
  assert.equal(await page.locator("#tooltip").isVisible(), false);
  await page.mouse.click(left.screenX, left.screenY, { button: "right" });
  assert.equal(await page.locator("#readMenu").isVisible(), true);
  assert.equal(await page.locator("#tooltip").isVisible(), false);
  await page
    .getByRole("menuitem", { name: "Read details", exact: true })
    .click();
  assert.match(await page.textContent("#tooltip"), /MAPQ: 60/);
  await page.click("#clearSelection", { force: true });
  const p = await point("left");
  await page.mouse.click(p.screenX, p.screenY, { button: "right" });
  await page
    .getByRole("menuitem", { name: "Enable hover details", exact: true })
    .click();
  await page.mouse.move(p.screenX + 2, p.screenY);
  await page.waitForFunction(() => !document.querySelector("#tooltip").hidden);
  await page.mouse.click(p.screenX, p.screenY, { button: "right" });
  await page
    .getByRole("menuitem", { name: "Disable hover details", exact: true })
    .click();
  await page.mouse.move(p.screenX + 3, p.screenY);
  assert.equal(await page.locator("#tooltip").isVisible(), false);
  // Shift-click two pairs sharing the same lane selects those two pairs,
  // without accidentally selecting the overlapping pair on the next lane.
  const clickPair = async (name, shift = false) => {
    const p = await point(name);
    if (shift) await page.keyboard.down("Shift");
    await page.mouse.click(p.screenX, p.screenY);
    if (shift) await page.keyboard.up("Shift");
  };
  await clickPair("left");
  await clickPair("right", true);
  assert.equal(await page.textContent("#selected"), "2 pairs selected");
  await page.fill("#hapName", "Manual H1");
  await page.click("#merge");
  assert.match(await page.textContent("#counts"), /1 manual haplotypes/);
  await page.click("#toggleAnalysis");
  assert.equal(await page.locator("#inspector").isVisible(), true);
  assert.match(
    await page.textContent("#inspector"),
    /2 included \/ 2 assigned pairs/,
  );
  assert.match(
    await page.textContent("#inspector"),
    /Site VAF 1\/1 \(100.0%\)/,
  );
  assert.match(await page.textContent("#treeSummary"), /1 manual haplotypes/);
  const narrower = (await canvas.boundingBox()).width;
  assert.ok(narrower < beforeWidth - 200);
  await page.click("#toggleAnalysis");
  assert.ok((await canvas.boundingBox()).width > narrower + 200);
  await page.evaluate(() => new Promise(requestAnimationFrame));
  const insertion = await page.evaluate(() =>
    window.__rects.find((r) => r.color === "#7618dc" && r.h >= 14),
  );
  assert.ok(insertion);
  const insertionBox = await canvas.boundingBox();
  await page.mouse.move(
    insertionBox.x + insertion.x + insertion.w / 2,
    insertionBox.y + insertion.y + insertion.h / 2,
  );
  await page.waitForFunction(() => !document.querySelector("#tooltip").hidden);
  assert.match(await page.textContent("#tooltip"), /ACGT/);
  await page.mouse.click(
    insertionBox.x + insertion.x + insertion.w / 2,
    insertionBox.y + insertion.y + insertion.h / 2,
  );
  assert.match(await page.textContent("#tooltip"), /ACGT/);
  assert.equal(await page.textContent("#selected"), "0 pairs selected");
  await page.click("#viewRange");
  // Group markers and raw markers have the same original CIGAR columns.
  await page.evaluate(() => new Promise(requestAnimationFrame));
  await page.waitForFunction(
    () =>
      Math.abs(
        document.querySelector("#pileup").width / devicePixelRatio -
          document.querySelector("#pileup").clientWidth,
      ) < 1,
  );
  const pos = await page.evaluate(() => {
    const c = document.querySelector("#pileup"),
      left = +c.dataset.leftMargin,
      x =
        left +
        ((1050 - +c.dataset.viewStart) /
          (+c.dataset.viewEnd - +c.dataset.viewStart)) *
          (c.clientWidth - left - 12);
    return { x, rects: window.__rects };
  });
  assert.ok(
    pos.rects.some((r) => r.color === "#0000ff" && Math.abs(r.x - pos.x) < 0.1),
  );
  await page.click("#toggleAnalysis");
  await page
    .getByRole("button", { name: "Expand read pairs", exact: true })
    .click();
  await page.click("#toggleAnalysis");
  const a = await info("left"),
    b = await info("right");
  assert.equal(a.y, b.y);
  await clickPair("left");
  await clickPair("right", true);
  await page.click("#merge");
  assert.equal(await page.textContent("#selected"), "0 pairs selected");
  await page.locator("#accessible").evaluate((e) => (e.open = true));
  assert.equal(
    await page
      .locator("#accessibleRows button")
      .filter({ hasText: /^left(?: ✓)?$/ })
      .count(),
    0,
  );
  await page.locator("#accessible").evaluate((e) => (e.open = false));
  await page.click("#toggleAnalysis");
  await page
    .getByRole("button", { name: "Expand read pairs", exact: true })
    .click();
  await page.click("#toggleAnalysis");

  await page.fill("#viewRange", "synthetic:1041-1060");
  await page.evaluate(() => {
    window.__texts = [];
    document.querySelector("#readScroll").style.width = "calc(100% - 17px)";
  });
  await page.click("#go");
  await page.evaluate(() => new Promise(requestAnimationFrame));
  const textAlignment = await page.evaluate(() => {
    const c = document.querySelector("#pileup"),
      left = +c.dataset.leftMargin,
      expected =
        left +
        ((1050 - +c.dataset.viewStart + 0.5) /
          (+c.dataset.viewEnd - +c.dataset.viewStart)) *
          (c.clientWidth - left - 12);
    const ts = window.__texts || [];
    return {
      ref: ts.some(
        (t) =>
          t.id === "referenceTrack" &&
          t.text === "A" &&
          t.align === "center" &&
          Math.abs(t.x - expected) < 0.1,
      ),
      read: ts.some(
        (t) =>
          t.id === "pileup" &&
          t.text === "C" &&
          t.align === "center" &&
          Math.abs(t.x - expected) < 0.1,
      ),
    };
  });
  assert.ok(
    textAlignment.ref && textAlignment.read,
    JSON.stringify(textAlignment),
  );
  await page.evaluate(
    () => (document.querySelector("#readScroll").style.width = ""),
  );
  await page.click("#fit");
  await page.evaluate(() => new Promise(requestAnimationFrame));
  // Hide a selected pair with no destructive removal; undo and restore retain it.
  await clickPair("overlap");
  await page.click("#hide");
  assert.match(await page.textContent("#counts"), /1 hidden pairs/);
  await page.click("#undo");
  assert.match(await page.textContent("#counts"), /0 hidden pairs/);
  await page.check("#includeDuplicates");
  assert.match(await page.textContent("#rowCount"), /5 visible pairs/);
  await page.uncheck("#includeDuplicates");
  assert.match(await page.textContent("#rowCount"), /4 visible pairs/);
  await page.selectOption("#readSize", "24");
  await page.waitForFunction(
    () => document.querySelector("#pileup").dataset.readHeight === "24",
  );
  assert.ok(
    await page.evaluate(() =>
      window.__rects.some((r) => r.h === 24 && r.color === "#c8c8c8"),
    ),
  );
  // Each drawer and the top panel can be collapsed independently.
  await page.click("#toggleInputs");
  assert.equal(await page.locator("#bam").isVisible(), true);
  await page.click("#toggleInputs");
  const compactHeight = (await canvas.boundingBox()).height;
  await page.click("#toggleTop");
  assert.equal(await page.locator(".workspace-head").isVisible(), true);
  await page.waitForFunction(
    () =>
      Math.abs(
        document.querySelector("#pileup").clientHeight -
          document.querySelector("#readScroll").clientHeight,
      ) < 2,
  );
  assert.ok((await canvas.boundingBox()).height < compactHeight);
  await page.click("#toggleTop");
  await page.screenshot({
    path: "/tmp/phase-packed-desktop.png",
    fullPage: true,
  });
  const exported = page.waitForEvent("download");
  await page.click("#toggleInputs");
  await page.click("#save");
  assert.equal(
    JSON.parse(await readFile(await (await exported).path(), "utf8")).groups
      .length,
    1,
  );
  await page.click("#toggleInputs");
  // Exercise the existing deep-coverage demo and a large somatic mixture with
  // no implicit merging, as well as reference/missing-SEQ binary regressions.
  for (const [name, index, ref] of [
    ["record-compatibility", "csi", "2bit"],
    ["somatic-ten", "csi", "fa"],
  ]) {
    await page.click("#reset");
    await page.setInputFiles("#bam", resolve(`tests/fixtures/${name}.bam`));
    await page.setInputFiles(
      "#index",
      resolve(`tests/fixtures/${name}.${index}`),
    );
    await page.setInputFiles(
      "#reference",
      resolve(`tests/fixtures/${name}.${ref}`),
    );
    if (ref === "fa")
      await page.setInputFiles(
        "#fai",
        resolve(`tests/fixtures/${name}.fa.fai`),
      );
    await page.click("#open");
    await page.waitForFunction(() =>
      document
        .querySelector("#status")
        .textContent.includes("Opened indexed BAM"),
    );
    await page.fill("#region", "synthetic:1001-1800");
    await page.click("#load");
    await page.waitForFunction(() =>
      document.querySelector("#status").textContent.includes("Loaded "),
    );
    assert.match(await page.textContent("#counts"), /0 manual haplotypes/);
    if (name === "somatic-ten")
      assert.match(
        await page.textContent("#counts"),
        /179 pairs\/fragments · 353 records/,
      );
  }
  await page.click("#zoomIn");
  const oldStart = +(await canvas.getAttribute("data-view-start")),
    box = await canvas.boundingBox();
  await page.mouse.move(box.x + box.width * 0.7, box.y + 80);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.7 - 70, box.y + 80, { steps: 5 });
  await page.mouse.up();
  assert.ok(+(await canvas.getAttribute("data-view-start")) > oldStart);
  await page
    .locator("#readScroll")
    .evaluate((e) => (e.scrollTop = e.scrollHeight));
  await page.waitForFunction(
    () => +document.querySelector("#pileup").dataset.firstRow > 0,
  );
  await page.click("#reset");
  assert.equal(await page.textContent("#tooltip"), "");
  assert.equal(await page.locator("#readMenu").isVisible(), false);
  assert.equal(await page.inputValue("#bam"), "");
  await page.fill("#gene", "TP53");
  await page.click("#find");
  const gene = await page.inputValue("#region");
  await page.selectOption("#exon", "0");
  assert.notEqual(await page.inputValue("#region"), gene);
  await page.setViewportSize({ width: 720, height: 1000 });
  await page.click("#toggleInputs");
  await page.screenshot({
    path: "/tmp/phase-packed-mobile.png",
    fullPage: true,
  });
  assert.ok(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth + 1,
    ),
  );
  assert.deepEqual(errors, []);
  assert.deepEqual(consoleErrors, []);
  assert.deepEqual(requests, []);
  console.log(
    "PASS: packed shared lanes with pair-specific hit testing, Shift selection, manual merges, consensus geometry, VAF/hierarchy, duplicate filtering, large glyphs, hover off and right-click details, optional hover, independent drawers/top compression, viewport sizing, BAM/BAI/CSI/2bit/FASTA, omitted SEQ, somatic depth, panning, scrolling, reset, responsive layout, no network/CSP errors.",
  );
} finally {
  await browser.close();
}
