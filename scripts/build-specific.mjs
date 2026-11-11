import { build } from "esbuild";
import { readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import path from "node:path";
const worker = await build({
  entryPoints: ["specific/worker.js"],
  bundle: true,
  write: false,
  format: "iife",
  platform: "browser",
  minify: true,
  plugins: [
    {
      name: "javascript-inflater",
      setup(b) {
        b.onResolve({ filter: /^@gmod\/bgzf-filehandle$/ }, () => ({
          path: path.resolve("src/bgzf.js"),
        }));
      },
    },
  ],
});
const app = await build({
  entryPoints: ["specific/app.js"],
  bundle: true,
  write: false,
  format: "iife",
  platform: "browser",
  minify: true,
  define: { WORKER_SOURCE: JSON.stringify(worker.outputFiles[0].text) },
});
const js = app.outputFiles[0].text,
  css = await readFile("specific/styles.css", "utf8");
if (/WebAssembly\s*\.|\.wasm["']|instantiateStreaming/.test(js))
  throw Error("WASM is not allowed in this release");
if (/<\/script/i.test(js)) throw Error("Unexpected script terminator");
const hash = (s) => createHash("sha256").update(s).digest("base64");
const policy = `default-src 'none'; script-src 'sha256-${hash(js)}'; script-src-attr 'none'; style-src 'sha256-${hash(css)}'; style-src-attr 'none'; worker-src blob:; connect-src 'none'; img-src 'none'; font-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'; frame-src 'none'; media-src 'none'; manifest-src 'none'`;
const html = (await readFile("specific/template.html", "utf8"))
  .replace(
    "<!--CSP-->",
    `<meta http-equiv="Content-Security-Policy" content="${policy}">`,
  )
  .replace("<!--STYLE-->", () => `<style>${css}</style>`)
  .replace("<!--SCRIPT-->", () => `<script>${js}</script>`);
await writeFile("specific/index.html", html);
await writeFile(
  "specific/SHA256SUMS",
  `${createHash("sha256").update(html).digest("hex")}  index.html\n`,
);
const notices = [];
for (const name of [
  "@gmod/bam",
  "@gmod/shared-read-cache",
  "@jbrowse/quick-lru",
  "generic-filehandle2",
  "pako",
]) {
  const pkg = JSON.parse(
    await readFile(`node_modules/${name}/package.json`, "utf8"),
  );
  let license = "";
  for (const file of ["LICENSE", "LICENSE.txt", "LICENSE.md"]) {
    try {
      license = await readFile(`node_modules/${name}/${file}`, "utf8");
      break;
    } catch {}
  }
  notices.push(`${name} ${pkg.version} (${pkg.license})\n${license}`);
}
await writeFile("specific/THIRD_PARTY_NOTICES.txt", notices.join("\n\n"));
console.log(
  `Built specific/index.html (${(Buffer.byteLength(html) / 1024 / 1024).toFixed(1)} MiB), standalone offline HTML with CSP hashes`,
);
