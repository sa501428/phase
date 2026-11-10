import { TwoBitFile } from "./twobit.js";
export class LocalFile {
  constructor(file, budget = { bytes: 0 }) {
    this.file = file;
    this.budget = budget;
  }
  async read(length, position = 0) {
    if (
      !Number.isSafeInteger(length) ||
      !Number.isSafeInteger(position) ||
      position < 0 ||
      length < 0 ||
      length > 64 * 1024 * 1024
    )
      throw Error("Invalid or excessive file range");
    this.budget.bytes += Math.min(
      length,
      Math.max(0, this.file.size - position),
    );
    if (this.budget.bytes > 256 * 1024 * 1024)
      throw Error("File I/O budget exceeded; narrow the interval");
    return new Uint8Array(
      await this.file.slice(position, position + length).arrayBuffer(),
    );
  }
  async readFile(options) {
    if (this.file.size > 64 * 1024 * 1024)
      throw Error("Index/text file exceeds 64 MiB limit");
    const a = await this.read(this.file.size);
    return typeof options === "string" || options?.encoding
      ? new TextDecoder().decode(a)
      : a;
  }
  async stat() {
    return { size: this.file.size };
  }
  async close() {}
}
export async function referenceReader(file, fai) {
  const handle = new LocalFile(file);
  if (/\.2bit$/i.test(file.name)) {
    const t = new TwoBitFile({ filehandle: handle });
    return {
      size: (n) => t.getSequenceSize(n),
      sequence: async (n, s, e) =>
        (await t.getSequence(n, s, e))?.toUpperCase(),
      reset: () => (handle.budget.bytes = 0),
    };
  }
  if (!fai) throw Error("Select the matching .fai for an uncompressed FASTA");
  if (fai.size > 8e6) throw Error("FAI too large");
  const entries = new Map();
  for (const line of (await fai.text()).trim().split(/\r?\n/)) {
    const [n, ...v] = line.split("\t"),
      [length, offset, bases, width] = v.map(Number);
    if (
      !n ||
      ![length, offset, bases, width].every(Number.isSafeInteger) ||
      length < 0 ||
      offset < 0 ||
      bases < 1 ||
      width < bases
    )
      throw Error("Invalid FAI");
    entries.set(n, { length, offset, bases, width });
  }
  return {
    reset: () => (handle.budget.bytes = 0),
    size: async (n) => entries.get(n)?.length,
    sequence: async (n, s, e) => {
      const a = entries.get(n);
      if (!a || s < 0 || e > a.length)
        throw Error("Reference interval unavailable");
      const from = a.offset + Math.floor(s / a.bases) * a.width + (s % a.bases),
        to =
          a.offset +
          Math.floor((e - 1) / a.bases) * a.width +
          ((e - 1) % a.bases) +
          1;
      return new TextDecoder()
        .decode(await handle.read(to - from, from))
        .replace(/\s/g, "")
        .toUpperCase();
    },
  };
}
export async function resolveReference(reader, name, length) {
  const names = [name, name.startsWith("chr") ? name.slice(3) : "chr" + name];
  if (name === "MT" || name === "chrM") names.push("MT", "chrM");
  const found = [];
  for (const n of new Set(names)) {
    const size = await reader.size(n);
    if (size !== undefined) found.push([n, size]);
  }
  const exact = found.find((a) => a[0] === name),
    chosen = exact || (found.length === 1 && found[0]);
  if (!chosen) throw Error("Reference contig unavailable or ambiguous");
  if (chosen[1] !== length) throw Error("BAM/reference contig lengths differ");
  return chosen[0];
}
