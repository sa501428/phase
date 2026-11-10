import { chromium } from "playwright-core";
import { resolve } from "node:path";
import assert from "node:assert/strict";
import { demo } from "../v2/demo.js";
const browser = await chromium.launch({
  executablePath:
    process.env.CHROME_PATH ||
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  headless: true,
  args: ["--disable-background-networking"],
});
try {
  const page = await browser.newPage({
      viewport: { width: 1440, height: 1100 },
    }),
    errors = [],
    network = [],
    violations = [];
  page.on("pageerror", (e) => errors.push(e.message));
  page.on("request", (r) => {
    if (/^https?:/.test(r.url())) network.push(r.url());
  });
  page.on("console", (m) => {
    if (m.type() === "error") violations.push(m.text());
  });
  await page.goto("file://" + resolve("index.html"));
  assert.equal(await page.title(), "Phase · local clone viewer");
  assert.equal(
    await page.getAttribute("header a", "href"),
    "archive/index.html",
  );
  await page.click("#demo");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("104 input records retained"),
  );
  assert.match(await page.textContent("#blocks"), /Block 1 · 3 supported/);
  assert.match(await page.textContent("#blocks"), /Block 2 · 2 supported/);
  assert.equal(await page.locator(".path").count(), 5);
  await page.locator(".path").first().click();
  assert.match(
    await page.textContent("#inspection"),
    /12 exclusively assigned fragments \/ 36/,
  );
  await page.click("#expand");
  assert.match(
    await page.textContent("#read-count"),
    /12 fragments · 24 alignment records/,
  );
  await page.locator("#reads").evaluate((e) => (e.scrollTop = e.scrollHeight));
  await page.waitForFunction(
    () => +document.querySelector("#read-canvas").dataset.firstRow > 0,
  );
  await page.locator("#read-canvas").scrollIntoViewIfNeeded();
  const canvas = page.locator("#read-canvas"),
    box = await canvas.boundingBox(),
    before = +(await canvas.getAttribute("data-view-start"));
  await page.mouse.move(box.x + box.width / 2, box.y + 100);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width / 2 + 75, box.y + 100, { steps: 5 });
  await page.mouse.up();
  assert.ok(+(await canvas.getAttribute("data-view-start")) < before);
  await page.click("#allreads");
  assert.match(
    await page.textContent("#read-count"),
    /52 fragments · 104 alignment records/,
  );
  const d = demo(),
    r1 = d.reference[120],
    r2 = d.reference[160],
    a1 = r1 === "A" ? "T" : "A",
    a2 = r2 === "C" ? "G" : "C";
  await page.fill(
    "#variants",
    `synthetic:1121${r1}>${a1}\nsynthetic:1161${r2}>${a2}`,
  );
  await page.click("#query");
  await page.waitForFunction(() =>
    document
      .querySelector("#query-results")
      .textContent.includes("Both ALT: 24"),
  );
  assert.match(await page.textContent("#query-results"), /Cis supported/);
  await page.screenshot({ path: "/tmp/phase-v2-desktop.png", fullPage: true });
  await page.click("#close-reads");
  await page.uncheck("#demo-ref");
  await page.click("#demo");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#status")
        .textContent.includes("104 input records retained") &&
      !document.querySelector("#demo").disabled,
  );
  await page.click("#query");
  await page.waitForFunction(() =>
    document
      .querySelector("#query-results")
      .textContent.includes("REF unverified"),
  );
  assert.match(await page.textContent("#query-results"), /Both ALT: 24/);
  const dl = page.waitForEvent("download");
  await page.click("#gfa");
  assert.equal((await dl).suggestedFilename(), "phase-local-assembly.gfa");
  await page.selectOption("#demo-kind", "somatic");
  await page.click("#demo");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("353 input records retained"),
  );
  assert.match(
    await page.textContent("#blocks"),
    /10 supported local patterns/,
  );
  assert.match(await page.textContent("#blocks"), /5 shared\/ambiguous/);
  assert.equal(await page.locator(".path").count(), 10);
  await page.screenshot({ path: "/tmp/phase-v2-somatic.png", fullPage: true });
  // Exercise actual binary BAM indexes and both reference formats, including
  // legal omitted SEQ and malformed per-record display geometry.
  for (const [fixture, index, reference] of [
    ["synthetic", "bai", null],
    ["paired", "csi", "fa"],
    ["record-compatibility", "bai", "2bit"],
    ["structural", "bai", null],
    ["somatic-ten", "csi", null],
  ]) {
    await page.click("#reset");
    await page.setInputFiles("#bam", resolve(`tests/fixtures/${fixture}.bam`));
    await page.setInputFiles(
      "#index",
      resolve(`tests/fixtures/${fixture}.${index}`),
    );
    if (reference) {
      await page.setInputFiles(
        "#reference",
        resolve(`tests/fixtures/${fixture}.${reference}`),
      );
      if (reference === "fa")
        await page.setInputFiles(
          "#fai",
          resolve(`tests/fixtures/${fixture}.fa.fai`),
        );
    }
    await page.click("#open");
    await page.waitForFunction(() =>
      document.querySelector("#status").textContent.includes("BAM opened"),
    );
    await page.fill("#chr", "synthetic");
    await page.fill("#start", "1001");
    await page.fill("#end", "1800");
    await page.click("#run");
    await page.waitForFunction(() =>
      document
        .querySelector("#status")
        .textContent.includes("input records retained"),
    );
    if (fixture === "record-compatibility") {
      assert.match(
        await page.textContent("#status"),
        /50 input records retained · 44 quality-eligible/,
      );
      assert.match(await page.textContent("#graph-note"), /3 missing-sequence/);
    }
    if (fixture === "structural")
      assert.match(
        await page.textContent("#structural"),
        /translocation|inversion/i,
      );
    if (fixture === "somatic-ten") {
      assert.match(
        await page.textContent("#blocks"),
        /10 supported local patterns/,
      );
      assert.match(await page.textContent("#blocks"), /5 shared\/ambiguous/);
    }
    if (fixture === "paired") {
      await page.fill("#start", "1051");
      await page.fill("#end", "1110");
      await page.click("#run");
      await page.waitForFunction(() =>
        document
          .querySelector("#status")
          .textContent.includes(
            "24 input records retained · 12 quality-eligible",
          ),
      );
      assert.match(await page.textContent("#graph-note"), /12 outside-target/);
    }
    await page
      .locator("#segments")
      .evaluate((e) => (e.parentElement.open = true));
    await page.locator("#segments button").first().click();
    await page.click("#expand");
    await page.click("#allreads");
    assert.equal(await page.locator("#read-canvas").count(), 1);
  }
  await page.click("#reset");
  assert.equal(await page.inputValue("#bam"), "");
  assert.equal(await page.locator(".path").count(), 0);
  assert.equal(await page.locator("#read-canvas").count(), 0);
  await page
    .locator("#gene")
    .evaluate((e) => (e.closest("details").open = true));
  await page.fill("#gene", "TP53");
  await page.click("#find");
  assert.match(await page.inputValue("#chr"), /17/);
  assert.ok(
    +(await page.inputValue("#end")) > +(await page.inputValue("#start")),
  );
  const transcriptSpan =
    +(await page.inputValue("#end")) - +(await page.inputValue("#start"));
  await page.selectOption("#exon", "0");
  assert.ok(
    +(await page.inputValue("#end")) - +(await page.inputValue("#start")) <
      transcriptSpan,
  );
  await page.setViewportSize({ width: 720, height: 1000 });
  await page.screenshot({ path: "/tmp/phase-v2-mobile.png", fullPage: true });
  assert.ok(
    await page.evaluate(
      () => document.documentElement.scrollWidth <= innerWidth + 1,
    ),
  );
  assert.deepEqual(network, []);
  assert.deepEqual(errors, []);
  assert.deepEqual(violations, []);
  // The frozen version remains independently runnable at the requested URL.
  await page.goto("file://" + resolve("archive/index.html"));
  await page.click("#demo2");
  await page.waitForFunction(() =>
    document
      .querySelector("#status")
      .textContent.includes("eligible fragments analyzed"),
  );
  assert.match(await page.textContent("#blockChips"), /B1: 3 groups/);
  assert.deepEqual(network, []);
  assert.deepEqual(errors, []);
  console.log(
    "PASS: v2 file:// graph, multiple haplotypes, paired reads, scroll, pan, queries, GFA, indexed BAM/BAI/CSI, FASTA/2bit, missing SEQ, structural clues, reset, offline catalog, responsive layout, no network/CSP errors, and archived viewer.",
  );
} finally {
  await browser.close();
}
