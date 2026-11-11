import { chromium } from "playwright-core";
import { resolve } from "node:path";
import { mkdir, writeFile } from "node:fs/promises";
import assert from "node:assert/strict";
import { fixture } from "./fixtures.mjs";
const dir = "/tmp/phase-specific-tests";
await mkdir(dir, { recursive: true });
const seq = (start, length, changes) => {
  const s = Array(length).fill("A");
  for (const [p, b] of changes)
    if (p >= start && p < start + length) s[p - start] = b;
  return s.join("");
};
function pair(name, x, y, changes) {
  return [
    [x, y, 99],
    [y, x, 147],
  ].map(([start, mate, flags]) => ({
    name,
    start,
    mate,
    flags,
    mateRef: 0,
    mq: 60,
    seq: seq(start, 50, changes),
    cigar: "50M",
    qual: Array(50).fill(35),
  }));
}
const direct = [];
let n = 0;
for (const [count, changes] of [
  [
    10,
    [
      [1120, "C"],
      [1700, "G"],
    ],
  ],
  [6, [[1120, "C"]]],
  [4, [[1700, "G"]]],
  [8, []],
])
  for (let i = 0; i < count; i++) {
    direct.push(...pair(`direct-${n}`, 1100 + n, 1680 + n, changes));
    n++;
  }
// Keep the target inside each end even as distinct start coordinates vary.
for (let i = 0; i < direct.length; i += 2) {
  const index = i / 2;
  const changes =
    index < 10
      ? [
          [1120, "C"],
          [1700, "G"],
        ]
      : index < 16
        ? [[1120, "C"]]
        : index < 20
          ? [[1700, "G"]]
          : [];
  const replacement = pair(direct[i].name, 1090 + index, 1670 + index, changes);
  direct[i] = replacement[0];
  direct[i + 1] = replacement[1];
}
const bridge = [];
for (const [block, x, y, changes] of [
  [
    0,
    1100,
    1180,
    [
      [1120, "C"],
      [1200, "G"],
    ],
  ],
  [
    1,
    1180,
    1580,
    [
      [1200, "G"],
      [1600, "T"],
    ],
  ],
  [
    2,
    1580,
    1680,
    [
      [1600, "T"],
      [1700, "C"],
    ],
  ],
])
  for (let i = 0; i < 12; i++)
    bridge.push(
      ...pair(`block${block}-${i}`, x + i, y + i, i < 6 ? changes : []),
    );
for (const [name, reads] of [
  ["direct", direct],
  ["bridge", bridge],
]) {
  const data = fixture({
    reference: "A".repeat(800),
    region: { chr: "synthetic", start: 1000, end: 1800 },
    reads,
  });
  for (const ext of ["bam", "bai", "csi"])
    await writeFile(`${dir}/${name}.${ext}`, data[ext]);
}
const browser = await chromium.launch({
  executablePath:
    process.env.CHROME_PATH ||
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  headless: true,
  args: ["--disable-background-networking"],
});
try {
  const page = await browser.newPage({
      viewport: { width: 1500, height: 1100 },
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
  await page.goto("file://" + resolve("specific/index.html"));
  assert.equal(await page.inputValue("#minVaf"), "2");
  await page.screenshot({ path: `${dir}/initial.png`, fullPage: true });
  async function open(name, index = "bai") {
    await page.setInputFiles("#bam", `${dir}/${name}.bam`);
    await page.setInputFiles("#index", `${dir}/${name}.${index}`);
    await page.click("#open");
    await page.waitForFunction(() =>
      document.querySelector("#fileStatus").textContent.includes("contigs"),
    );
  }
  async function check(text) {
    await page.fill("#variants", text);
    await page.click("#parse");
    await page.waitForFunction(
      () => !document.querySelector("#parse").disabled,
    );
  }
  async function run() {
    await page.click("#run");
    await page.waitForFunction(
      () =>
        document.querySelector("#status").textContent.startsWith("Complete"),
      {},
      { timeout: 60000 },
    );
  }
  await open("direct");
  await page.uncheck("#helpers");
  await check("synthetic:1121A>C\nsynthetic\t1701\tA\tG");
  await run();
  await page.fill("#minVaf", "100");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#runSettings")
        .textContent.includes("minimum VAF 100%") &&
      !document.querySelector("#run").disabled,
  );
  assert.equal(await page.locator("#pairRows tr[data-pair]").count(), 0);
  assert.match(
    await page.locator("#status").innerText(),
    /2 variants suppressed/,
  );
  assert.equal(await page.locator("#evidence").isVisible(), false);
  assert.equal(await page.locator("#evidenceRows tr").count(), 0);
  await page.fill("#minVaf", "0");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#runSettings")
        .textContent.includes("minimum VAF 0%") &&
      !document.querySelector("#run").disabled,
  );
  assert.equal(await page.locator("#variantA option").count(), 2);
  await page.fill("#minVaf", "2");
  await page.waitForFunction(
    () =>
      document
        .querySelector("#runSettings")
        .textContent.includes("minimum VAF 2%") &&
      !document.querySelector("#run").disabled,
  );
  let row = await page.locator("#pairRows tr").first().innerText();
  assert.match(row, /50\.0%/);
  assert.match(row, /10 \/ 20/);
  assert.match(row, /Direct read pairs/);
  assert.equal(await page.locator("#evidenceRows tr").count(), 28);
  assert.match(await page.locator("#evidenceRows").innerText(), /R1 synthetic/);
  assert.match(await page.locator("#evidenceRows").innerText(), /R2 synthetic/);
  assert.match(await page.locator("#pairDetail").innerText(), /neither 8/);
  await page.screenshot({ path: `${dir}/direct.png`, fullPage: true });
  // Invalid annotations block analysis; correcting the list clears old results.
  await check("synthetic:1121A>C\np.R1G");
  assert.equal(await page.locator("#run").isDisabled(), true);
  assert.equal(await page.locator("#results").isVisible(), false);
  assert.match(await page.locator("#parseResults").innerText(), /Protein-only/);
  // Real indexed file queries find a helper chain without loading a reference.
  await open("bridge", "csi");
  await page.check("#helpers");
  await check("synthetic:1121A>C\nsynthetic:1701A>C");
  await run();
  assert.match(
    await page.locator("#status").innerText(),
    /2 qualifying helper SNPs/,
  );
  assert.match(
    await page.locator("#pairDetail").innerText(),
    /No directly measured joint REF\/ALT observations/,
  );
  assert.match(
    await page.locator("#pathDetail").innerText(),
    /Helper-assisted phase link/,
  );
  assert.equal(await page.locator(".path-node.helper").count(), 2);
  assert.match(
    await page.locator("#pathDetail").innerText(),
    /12\/12 \(100\.0%\)/,
  );
  assert.equal(await page.locator("#helperRows tr.used-helper").count(), 2);
  await page.selectOption("#pairMode", "all");
  assert.match(
    await page.locator("#pairRows").innerText(),
    /Helper link available/,
  );
  await page.locator(".path-edges button").first().click();
  await page.waitForFunction(() =>
    document.querySelector("#evidenceRows").textContent.includes("block0"),
  );
  assert.equal(await page.locator("#evidenceRows tr").count(), 12);
  assert.match(await page.locator("#pairDetail").innerText(), /100\.0%/);
  await page.screenshot({ path: `${dir}/helpers.png`, fullPage: true });
  // Coding annotation uses the embedded offline RefSeq catalog, not a network API.
  await check("NM_000546.6:c.743G>A");
  assert.match(await page.locator("#parseResults").innerText(), /chr17:/);
  assert.equal(await page.locator("#run").isDisabled(), false);
  // Larger lists are paginated rather than cut at the historical 32-site cap.
  await check(
    Array.from({ length: 80 }, (_, i) => `synthetic ${1020 + i} A C`).join(
      "\n",
    ),
  );
  assert.match(await page.locator("#parseSummary").innerText(), /80 unique/);
  assert.equal(await page.locator(".parse-row").count(), 50);
  await page.locator("#parsePager button").last().click();
  assert.equal(await page.locator(".parse-row").count(), 30);
  // Clear is destructive only to app memory; inputs and visible evidence vanish.
  await page.click("#clear");
  assert.equal(await page.inputValue("#variants"), "");
  assert.equal(await page.locator("#results").isVisible(), false);
  assert.equal(await page.locator("#parseResults").innerText(), "");
  assert.equal(await page.locator("#bam").evaluate((e) => e.files.length), 0);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.screenshot({ path: `${dir}/mobile.png`, fullPage: true });
  assert.equal(
    await page.evaluate(
      () => document.documentElement.scrollWidth > window.innerWidth,
    ),
    false,
  );
  assert.deepEqual(errors, []);
  assert.deepEqual(requests, []);
  assert.deepEqual(consoleErrors, []);
  assert.equal(
    await page.evaluate(() => localStorage.length + sessionStorage.length),
    0,
  );
  console.log(
    "Specific browser checks passed: paired-end frequencies, BAM/BAI/CSI, helpers, provenance, coding HGVS, pagination, reset, mobile layout and zero network/storage.",
  );
} finally {
  await browser.close();
}
