// Both UCSC byte orders; version 0 indexes. Ranges are zero-based half-open.
export class TwoBitFile {
  constructor({ filehandle }) {
    this.f = filehandle;
    this.index = null;
    this.records = new Map();
  }
  async view(n, p) {
    const b = await this.f.read(n, p);
    if (b.length !== n) throw Error("Truncated two-bit reference");
    return new DataView(b.buffer, b.byteOffset, b.byteLength);
  }
  async init() {
    if (this.index) return;
    const h = await this.view(16, 0);
    this.le = h.getUint32(0, true) === 0x1a412743;
    if (h.getUint32(0, this.le) !== 0x1a412743)
      throw Error("Invalid two-bit signature");
    if (h.getUint32(4, this.le) !== 0)
      throw Error("Only version 0 two-bit references are supported");
    const n = h.getUint32(8, this.le);
    if (n > 100000) throw Error("Too many reference contigs");
    const index = new Map();
    let p = 16;
    for (let i = 0; i < n; i++) {
      const size = (await this.view(1, p)).getUint8(0);
      if (!size) throw Error("Empty reference contig");
      const bytes = await this.f.read(size, p + 1),
        name = new TextDecoder().decode(bytes),
        offset = (await this.view(4, p + 1 + size)).getUint32(0, this.le);
      index.set(name, offset);
      p += size + 5;
    }
    this.index = index;
  }
  async record(name) {
    await this.init();
    if (!this.index.has(name)) return undefined;
    if (this.records.has(name)) return this.records.get(name);
    let p = this.index.get(name);
    const h = await this.view(8, p),
      length = h.getUint32(0, this.le),
      count = h.getUint32(4, this.le);
    p += 8;
    const blocks = async (count) => {
      if (count > 1000000) throw Error("Two-bit block budget exceeded");
      const v = await this.view(count * 8, p);
      p += count * 8;
      const out = [];
      for (let i = 0; i < count; i++) {
        const start = v.getUint32(i * 4, this.le),
          size = v.getUint32((i + count) * 4, this.le);
        if (start + size > length) throw Error("Invalid reference block");
        out.push([start, start + size]);
      }
      return out;
    };
    const unknown = await blocks(count),
      m = (await this.view(4, p)).getUint32(0, this.le);
    p += 4;
    await blocks(m);
    p += 4;
    const r = { length, unknown, offset: p };
    this.records.set(name, r);
    return r;
  }
  async getSequenceSize(name) {
    return (await this.record(name))?.length;
  }
  async getSequence(name, start, end) {
    const r = await this.record(name);
    if (!r) return undefined;
    if (start < 0 || end > r.length || end < start)
      throw Error("Invalid reference interval");
    const first = Math.floor(start / 4),
      last = Math.ceil(end / 4),
      bytes = await this.f.read(last - first, r.offset + first);
    if (bytes.length !== last - first)
      throw Error("Truncated reference sequence");
    const out = Array(end - start);
    for (let p = start; p < end; p++)
      out[p - start] = "TCAG"[
        (bytes[Math.floor(p / 4) - first] >> ((3 - (p % 4)) * 2)) & 3
      ];
    for (const [s, e] of r.unknown)
      for (let p = Math.max(s, start); p < Math.min(e, end); p++)
        out[p - start] = "N";
    return out.join("");
  }
}
