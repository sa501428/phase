import { chromium } from "playwright-core";
import { resolve } from "node:path";
import assert from "node:assert/strict";
import { readFile, writeFile } from "node:fs/promises";
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
    network = [],
    consoleErrors = [];
  page.on("pageerror", (e) => errors.push(e.message));
  page.on("request", (r) => {
    if (/^https?:/.test(r.url())) network.push(r.url());
  });
  page.on("console", (m) => {
    if (m.type() === "error") consoleErrors.push(m.text());
  });
  await page.addInitScript(() => {
    const original = CanvasRenderingContext2D.prototype.fillRect;
    CanvasRenderingContext2D.prototype.fillRect = function (x, y, w, h) {
      if (this.canvas.id === "pileup") {
        if (x === 0 && y === 0 && h === 570) window.__rects = [];
        (window.__rects ??= []).push({ x, y, w, h, color: this.fillStyle });
      }
      return original.call(this, x, y, w, h);
    };
  });
  await page.goto("file://" + resolve("index.html"));
  assert.equal(await page.title(), "Phase · manual haplotype workspace");
  assert.equal(await page.locator('a[href="archive2/index.html"]').count(), 1);
  await page
    .locator("#demo")
    .evaluate((e) => (e.closest("details").open = true));
  await page.click("#demo");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("Loaded 52 read pairs"),
  );
  assert.match(await page.textContent("#counts"), /0 manual haplotypes/);
  assert.equal(await page.locator(".group-chip").count(), 0);
  assert.equal(
    await page.locator("#pileup").getAttribute("data-row-count"),
    "52",
  );
  await page.selectOption("#sort", "name");
  const canvas = page.locator("#pileup");
  await canvas.scrollIntoViewIfNeeded();
  const rowClick = async (index, shift = false, x = 65) => {
    await canvas.scrollIntoViewIfNeeded();
    await page.locator("#readScroll").evaluate((e, index) => {
      if (index * 42 < e.scrollTop || index * 42 + 35 > e.scrollTop + 570)
        e.scrollTop = index * 42;
    }, index);
    const scroll = await page
      .locator("#readScroll")
      .evaluate((e) => e.scrollTop);
    const box = await canvas.boundingBox();
    if (shift) await page.keyboard.down("Shift");
    await page.mouse.click(box.x + x, box.y + index * 42 - scroll + 20);
    if (shift) await page.keyboard.up("Shift");
  };
  // Real shift-click gestures cover multiple physical pairs, then merge them.
  await rowClick(16);
  await rowClick(27, true);
  assert.equal(await page.textContent("#selected"), "12 pairs selected");
  await page.fill("#hapName", "Founder reference");
  await page.click("#merge");
  assert.equal(await page.locator(".group-chip").count(), 1);
  assert.match(
    await page.textContent("#inspector"),
    /12 included \/ 12 assigned pairs/,
  );
  await rowClick(1);
  await rowClick(3, true);
  assert.equal(await page.textContent("#selected"), "3 pairs selected");
  await page.click("#hide");
  assert.match(await page.textContent("#counts"), /3 hidden pairs/);
  await page.click("#undo");
  assert.match(await page.textContent("#counts"), /0 hidden pairs/);
  await page.click("#redo");
  assert.match(await page.textContent("#counts"), /3 hidden pairs/);
  await page.click("#restore");
  assert.match(await page.textContent("#counts"), /0 hidden pairs/);
  await rowClick(17);
  await rowClick(28, true);
  await page.fill("#hapName", "Shared SNPs");
  await page.click("#merge");
  assert.equal(await page.locator(".group-chip").count(), 2);
  assert.match(
    await page.textContent("#inspector"),
    /Site VAF 24\/36 \(66.7%\)/,
  );
  assert.match(
    await page.textContent("#inspector"),
    /Within group 12\/12 \(100.0%\)/,
  );
  // A collapsed SNP uses the same original genomic column as its raw members.
  await page.locator("#readScroll").evaluate((e) => (e.scrollTop = 0));
  await page.waitForFunction(
    () => document.querySelector("#pileup").dataset.firstRow === "0",
  );
  const at = await page.evaluate(() => {
    const c = document.querySelector("#pileup"),
      left = c.clientWidth < 570 ? 152 : 210,
      x =
        left +
        ((1120 - +c.dataset.viewStart) /
          (+c.dataset.viewEnd - +c.dataset.viewStart)) *
          (c.clientWidth - left - 12);
    return { x, rects: window.__rects };
  });
  assert.ok(
    at.rects.some(
      (r) => r.color === "#27915e" && Math.abs(r.x - at.x) < 0.05 && r.y === 56,
    ),
  );
  // Conflicting reference/ALT groups cannot merge, even when explicitly selected.
  await rowClick(0);
  await rowClick(1, true);
  await page.click("#merge");
  assert.match(
    await page.textContent("#status"),
    /Cannot merge: contradictory/,
  );
  assert.equal(await page.locator(".group-chip").count(), 2);
  await page.click("#clearSelection");
  // Expand and verify MAPQ and legal soft clips survive in the original pair.
  await page.locator(".group-chip").first().click();
  await page
    .getByRole("button", { name: "Expand read pairs", exact: true })
    .click();
  await canvas.scrollIntoViewIfNeeded();
  const bb = await canvas.boundingBox();
  const x = await page.evaluate(() => {
    const c = document.querySelector("#pileup"),
      left = c.clientWidth < 570 ? 152 : 210;
    return (
      left +
      ((1080 - +c.dataset.viewStart) /
        (+c.dataset.viewEnd - +c.dataset.viewStart)) *
        (c.clientWidth - left - 12)
    );
  });
  await page.mouse.move(bb.x + x, bb.y + 42 + 15);
  await page.waitForFunction(() => !document.querySelector("#tooltip").hidden);
  assert.match(await page.textContent("#tooltip"), /MAPQ: 42/);
  assert.match(await page.textContent("#tooltip"), /4S220M3S/);
  assert.match(await page.textContent("#tooltip"), /R2 MAPQ 58/);
  await page.mouse.move(bb.x + bb.width - 10, bb.y + 5);
  await page
    .locator("#readScroll")
    .evaluate((e) => (e.scrollTop = e.scrollHeight));
  await page.waitForFunction(
    () => +document.querySelector("#pileup").dataset.firstRow > 0,
  );
  // Zoom and direct dragging use a shared reference/read axis.
  await page.click("#zoomIn");
  const before = +(await canvas.getAttribute("data-view-start"));
  await canvas.scrollIntoViewIfNeeded();
  const pan = await canvas.boundingBox();
  await page.mouse.move(pan.x + pan.width * 0.7, pan.y + 100);
  await page.mouse.down();
  await page.mouse.move(pan.x + pan.width * 0.7 - 60, pan.y + 100, {
    steps: 5,
  });
  await page.mouse.up();
  assert.ok(+(await canvas.getAttribute("data-view-start")) > before);
  await page.screenshot({ path: "/tmp/phase-v3-desktop.png", fullPage: true });
  // Save/restore edits and explicit sequence exports.
  const saved = page.waitForEvent("download");
  await page.click("#save");
  const d = await saved;
  const savedText = await readFile(await d.path(), "utf8");
  assert.equal(JSON.parse(savedText).groups.length, 2);
  await writeFile("/tmp/phase-v3-assignments.json", savedText);
  const fa = page.waitForEvent("download");
  await page.click("#fasta");
  assert.match(await readFile(await (await fa).path(), "utf8"), />Shared_SNPs/);
  await page.click("#undo");
  await page.setInputFiles("#import", "/tmp/phase-v3-assignments.json");
  await page.waitForFunction(
    () =>
      document.querySelector("#status").textContent.includes("assignments") ||
      document.querySelector("#status").textContent.includes("restored"),
  );
  assert.equal(await page.locator(".group-chip").count(), 2);
  // Actual BAM/BAI/CSI and 2bit/FASTA loading, missing SEQ, large somatic mixture.
  for (const [name, index, ref] of [
    ["paired", "bai", "2bit"],
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
    assert.equal(await page.locator(".group-chip").count(), 0);
    if (name === "paired") {
      assert.match(
        await page.textContent("#counts"),
        /12 pairs\/fragments · 24 records/,
      );
      await page.fill("#region", "synthetic:1051-1110");
      await page.click("#load");
      await page.waitForFunction(() =>
        document.querySelector("#status").textContent.includes("Loaded "),
      );
      assert.match(
        await page.textContent("#counts"),
        /12 pairs\/fragments · 24 records/,
      );
    }
    if (name === "record-compatibility") {
      assert.match(await page.textContent("#counts"), /50 records/);
      await page
        .locator("#readScroll")
        .evaluate((e) => (e.scrollTop = e.scrollHeight));
    }
    if (name === "somatic-ten")
      assert.match(
        await page.textContent("#counts"),
        /179 pairs\/fragments · 353 records/,
      );
  }
  await page.click("#reset");
  assert.equal(await page.inputValue("#bam"), "");
  assert.equal(await page.locator(".group-chip").count(), 0);
  await page.fill("#gene", "TP53");
  await page.click("#find");
  const full = await page.inputValue("#region");
  await page.selectOption("#exon", "0");
  assert.notEqual(await page.inputValue("#region"), full);
  await page
    .locator("#demo")
    .evaluate((e) => (e.closest("details").open = true));
  await page.selectOption("#scenario", "somatic");
  await page.click("#demo");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("Loaded 179 read pairs"),
  );
  assert.equal(await page.locator(".group-chip").count(), 0);
  await page.setViewportSize({ width: 720, height: 1000 });
  await page.screenshot({ path: "/tmp/phase-v3-mobile.png", fullPage: true });
  assert.ok(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth + 1,
    ),
  );
  assert.deepEqual(network, []);
  assert.deepEqual(errors, []);
  assert.deepEqual(consoleErrors, []);
  await page.goto("file://" + resolve("archive2/index.html"));
  await page.click("#demo");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("104 input records retained"),
  );
  assert.match(await page.textContent("#blocks"), /3 supported/);
  assert.deepEqual(errors, []);
  assert.deepEqual(network, []);
  console.log(
    "PASS: manual file:// viewer, real Shift-click merge/hide, undo/redo, conflict veto, paired MAPQ/soft clips, SNP geometry, VAF, hierarchy, scrolling/panning/zoom, JSON/FASTA, BAM/BAI/CSI/2bit/FASTA, missing SEQ, somatic demo, exon scope, responsive layout, privacy reset, no network/CSP errors, archive2.",
  );
} finally {
  await browser.close();
}
