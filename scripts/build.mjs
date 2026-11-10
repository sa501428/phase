import { build } from "esbuild";
import { readFile, writeFile, mkdir } from "node:fs/promises";
import { createHash } from "node:crypto";
import path from "node:path";
const plugins = [
  {
    name: "javascript-inflater",
    setup(b) {
      b.onResolve({ filter: /^@gmod\/bgzf-filehandle$/ }, () => ({
        path: path.resolve("src/bgzf.js"),
      }));
    },
  },
];
const worker = await build({
  entryPoints: ["v3/worker.js"],
  bundle: true,
  write: false,
  format: "iife",
  platform: "browser",
  minify: true,
  plugins,
  metafile: true,
});
const app = await build({
  entryPoints: ["v3/app.js"],
  bundle: true,
  write: false,
  format: "iife",
  platform: "browser",
  minify: true,
  define: { WORKER_SOURCE: JSON.stringify(worker.outputFiles[0].text) },
  metafile: true,
});
const js = app.outputFiles[0].text,
  css = await readFile("v3/styles.css", "utf8");
if (/WebAssembly\s*\.|\.wasm["']|instantiateStreaming/.test(js))
  throw Error("WASM in release bundle");
const hash = (s) => createHash("sha256").update(s).digest("base64");
const policy = `default-src 'none'; script-src 'sha256-${hash(js)}'; script-src-attr 'none'; style-src 'sha256-${hash(css)}'; style-src-attr 'none'; worker-src blob:; connect-src 'none'; img-src 'none'; font-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'; frame-src 'none'; media-src 'none'; manifest-src 'none'`;
const html = (await readFile("v3/index.html", "utf8"))
  .replace(
    "<!--CSP-->",
    `<meta http-equiv="Content-Security-Policy" content="${policy}">`,
  )
  .replace("<!--STYLE-->", () => `<style>${css}</style>`)
  .replace(
    "<!--SCRIPT-->",
    () => `<script>${js.replace(/<\/script/gi, "<\\/script")}</script>`,
  );
// Reject embedded script terminators rather than hashing a transformed script.
if (html.includes("<\\/script"))
  throw Error("Unexpected embedded script closing tag");
await mkdir("dist", { recursive: true });
await writeFile("index.html", html);
await writeFile(
  "dist/SHA256SUMS",
  `${createHash("sha256").update(html).digest("hex")}  ../index.html\n`,
);
await writeFile(
  "dist/build-inputs.json",
  JSON.stringify(
    {
      worker: Object.keys(worker.metafile.inputs),
      app: Object.keys(app.metafile.inputs),
    },
    null,
    2,
  ),
);
console.log(
  `Built index.html (${(html.length / 1024 / 1024).toFixed(1)} MiB), CSP hashes and checksum generated`,
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
await writeFile("dist/THIRD_PARTY_NOTICES.txt", notices.join("\n\n"));
