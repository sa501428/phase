import { chromium } from "playwright-core";
import { resolve } from "node:path";
import assert from "node:assert/strict";
const browser = await chromium.launch({
  executablePath:
    process.env.CHROME_PATH ||
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  headless: true,
  args: ["--disable-background-networking"],
});
try {
  const page = await browser.newPage();
  const errors = [],
    network = [];
  page.on("pageerror", (e) => errors.push(e.message));
  page.on("request", (r) => {
    if (/^https?:/.test(r.url())) network.push(r.url());
  });
  await page.addInitScript(() => {
    const original = CanvasRenderingContext2D.prototype.fillRect;
    CanvasRenderingContext2D.prototype.fillRect = function (x, y, w, h) {
      if (this.canvas.id === "pileup") {
        if (x === 0 && y === 0 && this.fillStyle === "#ffffff")
          window.__pileupRects = [];
        (window.__pileupRects ??= []).push({
          x,
          y,
          w,
          h,
          color: this.fillStyle,
        });
      }
      return original.call(this, x, y, w, h);
    };
  });
  await page.goto("file://" + resolve("index.html"));
  await page.click("#demo2");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#status")
        .textContent.includes("eligible fragments analyzed"),
    { timeout: 20000 },
  );
  assert.match(await page.textContent("#metrics"), /44/);
  assert.match(await page.textContent("#events"), /duplication/);
  assert.match(await page.textContent("#blockChips"), /B1: 3 groups/);
  assert.match(await page.textContent("#blockChips"), /B2: 1 groups/);
  assert.match(
    await page.textContent("#treeSummary"),
    /3 supported groups → 2 similarity clusters/,
  );
  await page.locator("#treeCanvas").hover({ position: { x: 180, y: 35 } });
  assert.match(
    await page.textContent("#groupInfo"),
    /12 assigned eligible fragments \/ 36/,
  );
  await page.locator("#treeCanvas").press("ArrowDown");
  assert.match(await page.textContent("#groupInfo"), /24\/36 \(66.7%\)/);
  await page.locator("#treeCut").focus();
  await page.keyboard.press("End");
  await page.locator("#treeCut").dispatchEvent("input");
  assert.match(
    await page.textContent("#treeSummary"),
    /→ 1 similarity clusters/,
  );
  await page.selectOption("#treeBlock", "2");
  assert.match(
    await page.textContent("#treeSummary"),
    /B2: 1 supported groups/,
  );
  await page.selectOption("#treeBlock", "1");
  await page.locator("#treeCut").evaluate((e) => (e.value = "25"));
  await page.locator("#treeCut").dispatchEvent("input");
  assert.equal(
    await page.locator("#evidenceDetails").getAttribute("open"),
    null,
  );

  assert.match(await page.textContent("#rowCount"), /4 group rows in view/);
  assert.equal(await page.locator("#collapse").isChecked(), true);
  assert.equal(await page.inputValue("#detail"), "overview");
  assert.equal(await page.inputValue("#layer"), "all");
  const markerPositions = () =>
    page.evaluate(() => {
      const canvas = document.querySelector("#pileup"),
        start = +canvas.dataset.viewStart,
        end = +canvas.dataset.viewEnd,
        width = Math.max(
          400,
          document.querySelector("#canvasWrap").clientWidth,
        ),
        scale = (width - 130) / (end - start);
      return (window.__pileupRects || [])
        .filter(
          (r) =>
            (r.h === 10 && ["#24865e", "#3075bc"].includes(r.color)) ||
            (r.w === 4 && r.h === 16 && r.color === "#9c70b3"),
        )
        .map((r) => ({
          color: r.color,
          pos: start + (r.x + (r.color === "#9c70b3" ? 2 : 0) - 115) / scale,
        }));
    });
  const collapsedMarkers = await markerPositions();
  assert.ok(
    collapsedMarkers.some(
      (m) => m.color === "#9c70b3" && Math.abs(m.pos - 1200) < 1e-6,
    ),
  );
  for (const m of collapsedMarkers)
    assert.ok(
      Math.abs(
        m.pos -
          (m.color === "#24865e" ? 1120 : m.color === "#3075bc" ? 1160 : 1200),
      ) < 1e-6,
    );

  await page.locator("#pileup").click({ position: { x: 40, y: 120 } });
  assert.equal(
    await page.locator('#readTable tr[aria-expanded="true"]').count(),
    1,
  );
  await page.click("#expandAll");
  assert.match(await page.textContent("#rowCount"), /of 48/);
  await page.locator("#pileup").focus();
  await page.keyboard.press("End");
  await page.waitForFunction(() =>
    document.querySelector("#rowCount").textContent.includes("29–48"),
  );
  assert.match(await page.textContent("#readTable"), /Read 44/);
  const expandedMarkers = await markerPositions();
  assert.ok(expandedMarkers.some((m) => m.color === "#9c70b3"));
  for (const m of expandedMarkers)
    assert.ok(
      Math.abs(
        m.pos -
          (m.color === "#24865e" ? 1120 : m.color === "#3075bc" ? 1160 : 1200),
      ) < 1e-6,
    );

  await page.click("#collapseAll");
  assert.match(await page.textContent("#rowCount"), /of 4/);
  // A drag moves genomic content, while a click still expands its family.
  await page.click("#zoomIn");
  await page.locator("#pileup").scrollIntoViewIfNeeded();
  const bounds = await page.locator("#pileup").boundingBox();
  const expandedBefore = await page
    .locator('#readTable tr[aria-expanded="true"]')
    .count();
  const rangeBefore = await page
    .locator("#pileup")
    .getAttribute("data-view-start");
  await page.mouse.move(bounds.x + 350, bounds.y + 120);
  await page.mouse.down();
  await page.mouse.move(bounds.x + 250, bounds.y + 120, { steps: 8 });
  await page.mouse.up();
  const rangeAfterLeft = await page
    .locator("#pileup")
    .getAttribute("data-view-start");
  assert.notEqual(rangeBefore, rangeAfterLeft, errors.join("\n"));
  assert.equal(
    await page.locator('#readTable tr[aria-expanded="true"]').count(),
    expandedBefore,
  );
  await page.mouse.move(bounds.x + 250, bounds.y + 120);
  await page.mouse.down();
  await page.mouse.move(bounds.x + 350, bounds.y + 120, { steps: 8 });
  await page.mouse.up();
  assert.equal(
    await page.locator("#pileup").getAttribute("data-view-start"),
    rangeBefore,
  );
  assert.equal(
    await page
      .locator("#pileup")
      .evaluate((c) => c.classList.contains("dragging")),
    false,
  );
  await page.click("#fit");
  await page.locator("#events button").first().click();
  await page.waitForFunction(() =>
    document.querySelector("#alignment").textContent.includes("Candidate"),
  );
  await page.screenshot({ path: "/tmp/phase-synthetic.png", fullPage: true });
  await page.click("#clear");
  assert.equal(await page.locator("#workspace").isVisible(), false);
  for (const [id, name] of [
    ["bam", "synthetic.bam"],
    ["index", "synthetic.bai"],
    ["reference", "synthetic.fa"],
    ["fai", "synthetic.fa.fai"],
  ])
    await page.setInputFiles("#" + id, resolve("tests/fixtures/" + name));
  await page.check("#confirm");
  await page.click("#open");
  await page.waitForFunction(() =>
    document.querySelector("#status").textContent.startsWith("Files open"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#run");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#status")
        .textContent.includes("eligible fragments analyzed"),
    { timeout: 20000 },
  );
  assert.match(await page.textContent("#metrics"), /44/);
  assert.match(await page.textContent("#events"), /duplication/);
  await page.click("#clear");
  for (const [id, name] of [
    ["bam", "synthetic.bam"],
    ["index", "synthetic.csi"],
    ["reference", "synthetic.2bit"],
  ])
    await page.setInputFiles("#" + id, resolve("tests/fixtures/" + name));
  await page.check("#confirm");
  await page.click("#open");
  await page.waitForFunction(() =>
    document.querySelector("#status").textContent.startsWith("Files open"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#run");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  assert.match(await page.textContent("#metrics"), /44/);
  await page.fill("#gene", "TP53");
  await page.click("#lookup");
  assert.match(await page.inputValue("#region"), /^chr17:/);
  await page.click("#clear");
  for (const [id, name] of [
    ["bam", "paired.bam"],
    ["index", "paired.bai"],
    ["reference", "paired.fa"],
    ["fai", "paired.fa.fai"],
  ])
    await page.setInputFiles("#" + id, resolve("tests/fixtures/" + name));
  await page.check("#confirm");
  await page.click("#open");
  await page.waitForFunction(() =>
    document.querySelector("#status").textContent.startsWith("Files open"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#run");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  assert.match(await page.textContent("#rowCount"), /1 group rows in view/);
  assert.match(
    await page.textContent("#rowCount"),
    /12 total fragments \/ 24 records retained/,
  );
  await page.click("#expandAll");
  assert.match(await page.textContent("#readTable"), /Pair/);
  assert.match(await page.textContent("#readTable"), /50M \/ 50M/);
  await page.fill("#region", "synthetic:1051-1110");
  await page.click("#run");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  assert.match(
    await page.textContent("#rowCount"),
    /12 total fragments \/ 24 records retained/,
  );
  await page.click("#clear");
  for (const [id, name] of [
    ["bam", "structural.bam"],
    ["index", "structural.bai"],
    ["reference", "structural.fa"],
    ["fai", "structural.fa.fai"],
  ])
    await page.setInputFiles("#" + id, resolve("tests/fixtures/" + name));
  await page.check("#confirm");
  await page.click("#open");
  await page.waitForFunction(() =>
    document.querySelector("#status").textContent.startsWith("Files open"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#run");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  assert.match(
    await page.textContent("#structural"),
    /Possible translocation · split alignment/,
  );
  assert.match(
    await page.textContent("#structural"),
    /Inversion-like split alignment/,
  );
  assert.match(
    await page.textContent("#structural"),
    /mate on another chromosome/,
  );
  assert.match(await page.textContent("#structural"), /chr8:5001/);
  await page.locator("#structural summary").first().click();
  await page.locator("#structural button").first().click();
  assert.match(
    await page.textContent("#structural"),
    /offline assembly catalog/,
  );
  await page.click("#clear");
  for (const [id, name] of [
    ["bam", "noisy.bam"],
    ["index", "noisy.bai"],
    ["reference", "noisy.fa"],
    ["fai", "noisy.fa.fai"],
  ])
    await page.setInputFiles("#" + id, resolve("tests/fixtures/" + name));
  await page.check("#confirm");
  await page.click("#open");
  await page.waitForFunction(() =>
    document.querySelector("#status").textContent.startsWith("Files open"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#run");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  const strictRows = Number(
    (await page.textContent("#rowCount")).match(/^(\d+) group rows/)[1],
  );
  assert.ok(strictRows > 1 && strictRows <= 200);
  assert.match(
    await page.textContent("#rowCount"),
    /200 total fragments \/ 200 records retained/,
  );
  await page.selectOption("#detail", "families");
  assert.match(await page.textContent("#rowCount"), /200 families in view/);
  await page.selectOption("#detail", "overview");
  assert.match(
    await page.textContent("#rowCount"),
    new RegExp("^" + strictRows + " group rows"),
  );
  await page.click("#expandAll");
  assert.match(
    await page.textContent("#rowCount"),
    new RegExp("of " + (200 + strictRows) + " ·"),
  );
  await page.locator("#pileup").focus();
  await page.keyboard.press("End");
  await page.waitForFunction(() =>
    document
      .querySelector("#rowCount")
      .textContent.match(/rows \d+–(\d+) of \1/),
  );
  assert.match(await page.textContent("#readTable"), /Read/);
  await page.click("#collapseAll");
  assert.match(
    await page.textContent("#rowCount"),
    new RegExp("of " + strictRows + " ·"),
  );
  assert.deepEqual(errors, []);
  assert.deepEqual(network, []);
  await page.reload();
  assert.equal(await page.locator("#workspace").isVisible(), false);
  console.log(
    "Standalone Chrome: synthetic demo, BAM/BAI+FASTA and BAM/CSI+two-bit queries, gene lookup, event alignment, reset/reload, no page errors or outbound requests passed",
  );
} finally {
  await browser.close();
}
