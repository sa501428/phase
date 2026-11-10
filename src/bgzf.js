// JavaScript-only replacement for the optional WASM-enabled GMOD inflater.
import { Inflate } from "pako";
export const MAX_BGZF_BLOCK_SIZE = 65536,
  POOL_MIN_DECOMPRESSED_BYTES = Infinity;
const MAX = 128 * 1024 * 1024;
export function scanBgzfBlocks(input, base = 0, last = Infinity) {
  const out = [];
  let p = 0,
    total = 0;
  while (p + 18 <= input.length) {
    if (
      input[p] !== 31 ||
      input[p + 1] !== 139 ||
      input[p + 2] !== 8 ||
      input[p + 3] !== 4
    )
      throw Error("Invalid BGZF header");
    const x = input[p + 10] | (input[p + 11] << 8);
    let q = p + 12,
      size = 0;
    while (q + 4 <= p + 12 + x) {
      const n = input[q + 2] | (input[q + 3] << 8);
      if (input[q] === 66 && input[q + 1] === 67 && n === 2)
        size = (input[q + 4] | (input[q + 5] << 8)) + 1;
      q += 4 + n;
    }
    if (size < 26) throw Error("Missing BGZF block size");
    if (p + size > input.length) break; // header range reads may end within a block
    const v = new DataView(input.buffer, input.byteOffset + p + size - 4, 4),
      n = v.getUint32(0, true);
    if (n > 65536 || (total += n) > MAX)
      throw Error("Decompression budget exceeded; narrow the interval");
    out.push({
      inputOffset: p,
      compressedSize: size,
      decompressedSize: n,
      filePosition: base + p,
    });
    if (base + p >= last) break;
    p += size;
  }
  return out;
}
function decode(input, b) {
  const inflater = new Inflate({ chunkSize: 16384 }),
    pieces = [];
  let total = 0;
  inflater.onData = (part) => {
    total += part.length;
    if (total > 65536) throw Error("BGZF expansion exceeds block budget");
    pieces.push(part);
  };
  inflater.push(
    input.subarray(b.inputOffset, b.inputOffset + b.compressedSize),
    false,
  );
  if (inflater.err || !inflater.ended || total !== b.decompressedSize)
    throw Error("Corrupt or truncated BGZF block");
  return join(pieces);
}
export async function unzip(input) {
  if (!input.length) return input;
  const bs = scanBgzfBlocks(input);
  return join(bs.map((b) => decode(input, b)));
}
function join(arr) {
  const n = arr.reduce((s, a) => s + a.length, 0);
  if (n > MAX) throw Error("Decompression budget exceeded");
  const out = new Uint8Array(n);
  let p = 0;
  for (const a of arr) {
    out.set(a, p);
    p += a.length;
  }
  return out;
}
export async function unzipChunkSegments(input, { minv, maxv }) {
  const bs = scanBgzfBlocks(input, minv.blockPosition, maxv.blockPosition),
    segments = [],
    cpositions = [],
    dpositions = [];
  let d = minv.dataPosition;
  for (let i = 0; i < bs.length; i++) {
    const b = bs[i],
      a = decode(input, b),
      last = b.filePosition >= maxv.blockPosition;
    cpositions.push(b.filePosition);
    dpositions.push(d);
    const start = i === 0 ? minv.dataPosition : 0,
      end = last ? Math.min(maxv.dataPosition + 1, a.length) : a.length;
    if (start < end) segments.push(a.subarray(start, end));
    d += Math.max(0, a.length - start);
    if (last) {
      cpositions.push(b.filePosition + b.compressedSize);
      dpositions.push(d);
      break;
    }
  }
  return { segments, cpositions, dpositions };
}
export async function unzipChunkSlice(input, chunk) {
  const r = await unzipChunkSegments(input, chunk);
  return { ...r, buffer: join(r.segments) };
}
