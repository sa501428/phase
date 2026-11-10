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
  assert.match(await page.textContent("#rowCount"), /4 families in view/);
  assert.equal(await page.locator("#collapse").isChecked(), true);
  assert.equal(await page.inputValue("#layer"), "all");
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
  await page.click("#collapseAll");
  assert.match(await page.textContent("#rowCount"), /of 4/);
  // A drag moves genomic content, while a click still expands its family.
  await page.click("#zoomIn");
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
  assert.notEqual(rangeBefore, rangeAfterLeft);
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
  assert.match(await page.textContent("#rowCount"), /1 families in view/);
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
