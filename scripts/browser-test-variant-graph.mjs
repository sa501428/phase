import { chromium } from "playwright-core";
import { mkdir, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import assert from "node:assert/strict";
import { fixture } from "./fixtures.mjs";
const dir = "/tmp/phase-variant-network";
await mkdir(dir, { recursive: true });
const reads = [];
for (let f = 0; f < 100; f++)
  for (const [start, mate, flags] of [
    [1000, 1100, 99],
    [1100, 1000, 147],
  ]) {
    const changes = [];
    if (f < 36) changes.push([1050, "C"]);
    if (f < 30) changes.push([1100, "G"]);
    if (f >= 30 && f < 35) changes.push([1150, "T"]);
    if (f === 35) changes.push([1200, "G"]);
    const seq = Array(350).fill("A");
    for (const [p, b] of changes)
      if (p >= start && p < start + 350) seq[p - start] = b;
    reads.push({
      id: reads.length,
      name: "network-" + f,
      start,
      end: start + 350,
      mate,
      mateRef: 0,
      flags,
      seq: seq.join(""),
      cigar: "350M",
      qual: Array(350).fill(35),
      mq: 60,
    });
  }
const bytes = fixture({
  reads,
  reference: "A".repeat(800),
  region: { chr: "synthetic", start: 1000, end: 1800 },
});
for (const ext of ["bam", "bai", "twobit"])
  await writeFile(
    `${dir}/network.${ext === "twobit" ? "2bit" : ext}`,
    bytes[ext],
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
      viewport: { width: 1600, height: 1150 },
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
  await page.goto("file://" + resolve("index.html"));
  await page.setInputFiles("#bam", `${dir}/network.bam`);
  await page.setInputFiles("#index", `${dir}/network.bai`);
  await page.setInputFiles("#reference", `${dir}/network.2bit`);
  await page.click("#open");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("Opened indexed BAM"),
  );
  await page.fill("#region", "synthetic:1001-1800");
  await page.click("#load");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "3",
  );
  assert.equal(await page.locator("#variantGraph").getAttribute("open"), "");
  assert.equal(await page.locator("#variantNetwork [data-node]").count(), 3);
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 2);
  const bubble = await page
    .locator('#variantNetwork [data-node="0"] circle')
    .boundingBox();
  assert.ok(
    bubble.width >= 35 && bubble.width < 60,
    "bubbles stay readable in screen pixels",
  );
  const weights = await page
    .locator("#variantNetwork line[data-weight]")
    .evaluateAll((es) =>
      es.map((e) => ({
        weight: +e.dataset.weight,
        width: +e.getAttribute("stroke-width"),
        shared: +e.dataset.shared,
      })),
    );
  const strong = weights.find((e) => e.shared === 30),
    weak = weights.find((e) => e.shared === 5);
  assert.equal(strong.weight, 0.3);
  assert.equal(weak.weight, 0.05);
  assert.ok(strong.width > 8 * weak.width);
  await page.locator('#variantNetwork [data-edge="1"]').focus();
  await page.keyboard.press("Enter");
  assert.match(await page.textContent("#graphDetail"), /Both ALT 5\/100/);
  await page
    .getByRole("button", { name: "Select 5 supporting pairs", exact: true })
    .click();
  assert.match(await page.textContent("#selected"), /^5 pairs selected$/);
  await page.locator('#variantNetwork [data-node="0"] circle').click();
  assert.match(await page.textContent("#graphDetail"), /ALT 36\/100/);
  await page
    .getByRole("button", { name: "Select 36 supporting pairs", exact: true })
    .click();
  assert.match(await page.textContent("#selected"), /^36 pairs selected$/);
  await page
    .getByRole("button", { name: "Show site in pileup", exact: true })
    .click();
  assert.match(await page.textContent("#siteLabel"), /1051/);
  await page.selectOption("#graphWeight", "linkage");
  assert.match(await page.textContent("#graphDefinition"), /left-only ALT/);
  await page.selectOption("#graphWeight", "frequency");
  const node = page.locator('#variantNetwork [data-node="2"]'),
    before = await node.getAttribute("transform"),
    box = await node.locator("circle").boundingBox();
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.down();
  await page.mouse.move(
    box.x + box.width / 2 + 60,
    box.y + box.height / 2 - 20,
    { steps: 8 },
  );
  await page.mouse.up();
  assert.notEqual(await node.getAttribute("transform"), before);
  const view = await page.getAttribute("#variantNetwork", "viewBox");
  await page.click("#graphIn");
  assert.notEqual(await page.getAttribute("#variantNetwork", "viewBox"), view);
  await page.click("#graphFit");
  await page.fill("#graphSupport", "10");
  await page.locator("#graphSupport").press("Tab");
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 1);
  await page.fill("#graphSupport", "1000");
  await page.locator("#graphSupport").press("Tab");
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 0);
  assert.equal(await page.locator("#variantNetwork [data-node]").count(), 3);
  assert.equal(await page.locator("#graphMessage").isVisible(), false);
  await page.fill("#graphSupport", "1");
  await page.locator("#graphSupport").press("Tab");
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 2);
  const openHeight = (await page.locator("#readScroll").boundingBox()).height;
  await page.locator("#variantGraph > summary").click();
  await page.waitForFunction(
    (h) => document.querySelector("#readScroll").clientHeight > h + 100,
    openHeight,
  );
  await page.locator("#variantGraph > summary").click();
  await page.click("#toggleInputs");
  await page.fill("#minVaf", "1");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "4",
  );
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 3);
  await page.fill("#minVaf", "5.1");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "2",
  );
  await page.fill("#minVaf", "5");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "3",
  );
  await page.click("#toggleInputs");
  await page.locator('#variantNetwork [data-node="2"] circle').click();
  await page
    .getByRole("button", { name: "Select 5 supporting pairs", exact: true })
    .click();
  await page.click("#hide");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "2",
  );
  await page.click("#undo");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "3",
  );
  await page.screenshot({ path: `${dir}/desktop.png`, fullPage: true });
  await page.setViewportSize({ width: 720, height: 1000 });
  await page.waitForTimeout(100);
  assert.ok((await page.locator("#variantNetwork").boundingBox()).width > 500);
  assert.ok(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth,
    ),
  );
  await page.screenshot({ path: `${dir}/mobile.png`, fullPage: true });
  await page.click("#reset");
  assert.equal(await page.locator("#variantNetwork [data-node]").count(), 0);
  assert.match(await page.textContent("#graphMessage"), /Load a region/);
  await page.setViewportSize({ width: 1600, height: 1150 });
  await page.fill("#minVaf", "2");
  await page.locator("details.demo > summary").click();
  await page.selectOption("#scenario", "somatic");
  await page.click("#demo");
  await page.waitForFunction(
    () => document.querySelector("#variantNetwork").dataset.nodes === "9",
  );
  await page.screenshot({ path: `${dir}/somatic.png`, fullPage: true });
  await page.locator("#variantGraph > summary").click();
  await page.click("#reset");
  assert.equal(await page.locator("#variantNetwork [data-node]").count(), 0);
  assert.equal(await page.locator("#variantNetwork [data-edge]").count(), 0);
  assert.deepEqual(errors, []);
  assert.deepEqual(network, []);
  assert.deepEqual(consoleErrors, []);
  console.log(
    "PASS: visible weighted variant graph, paired counts, large line contrast, node/edge inspection, supporting-pair selection, locus focus, dragging/zoom, support filter, collapse, live 1%/5% VAF, hide/undo, responsive layout, reset and offline CSP.",
  );
} finally {
  await browser.close();
}
