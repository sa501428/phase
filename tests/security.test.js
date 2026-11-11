import { test } from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { createHash } from "node:crypto";
test("release inline script/style match CSP hashes and contain no WASM", async () => {
  for (const release of [
    "index.html",
    "automatic/index.html",
    "assembly/index.html",
    "specific/index.html",
  ]) {
    const html = await readFile(release, "utf8");
    assert.match(html, /id="minVaf"[\s\S]{0,180}value="2"/);
    for (const tag of ["script", "style"]) {
      const content = html.match(
        new RegExp(`<${tag}>([\\s\\S]*?)<\\/${tag}>`),
      )[1];
      const hash = createHash("sha256").update(content).digest("base64");
      assert.ok(html.includes(`'sha256-${hash}'`));
    }
    assert.ok(html.includes("connect-src 'none'"));
    assert.doesNotMatch(html, /WebAssembly\s*\.|instantiateStreaming/);
  }
});
