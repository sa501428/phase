// Desktop IGV light-theme defaults, verified against igvteam/igv commit
// 811bf7042e0b05d124ff421ee7100eec4789d179: preferences.tab,
// alignment/AlignmentRenderer.java and alignment/BaseRenderer.java.
export const IGV_COLORS = Object.freeze({
  A: "#009600",
  C: "#0000ff",
  G: "#d17105",
  T: "#ff0000",
  N: "#808080",
});
export const IGV_INSERTION = "#7618dc",
  IGV_DELETION = "#000000",
  IGV_SKIP = "#96b8c8",
  IGV_ALIGNMENT = "#c8c8c8";
export function qualityAlpha(q) {
  if (q === undefined || q === 255 || q >= 20) return 1;
  const alpha =
    q < 5 ? 0.2 : Math.max(0.2, Math.min(1, 0.1 + (0.9 * (q - 5)) / 15));
  return Math.floor(alpha * 10 + 0.5) / 10;
}
export function drawIGVBase(c, base, x, y, width, height, quality) {
  const old = c.globalAlpha;
  c.globalAlpha *= qualityAlpha(quality);
  c.fillStyle = IGV_COLORS[base.toUpperCase()] || "#000000";
  const size = Math.min(Math.floor(width), height, 12);
  if (size >= 8) {
    c.font = `bold ${size}px Arial,sans-serif`;
    c.textAlign = "center";
    c.textBaseline = "middle";
    c.fillText(base, x + width / 2, y + height / 2);
    c.textAlign = "left";
    c.textBaseline = "alphabetic";
  } else c.fillRect(x, y, Math.max(1, width), height);
  c.globalAlpha = old;
}
export function drawIGVInsertion(c, x, y, height, length, ppb) {
  const wing = height > 10 ? 2 : height > 5 ? 1 : 0;
  c.fillStyle = IGV_INSERTION;
  c.font = "12px Arial,sans-serif";
  const label = String(length),
    maxWidth = length * ppb,
    textWidth = c.measureText(label).width + 4,
    show = length > 1 && height >= 8 && textWidth < maxWidth,
    width = length > 1 ? Math.max(2, Math.min(textWidth, maxWidth)) : 2,
    left = length > 1 ? x - 1 - Math.ceil(width / 2) : x;
  c.fillRect(left, y, width, height);
  if (wing) {
    c.fillRect(left - wing, y, width + 2 * wing, 2);
    c.fillRect(left - wing, y + height - 2, width + 2 * wing, 2);
  }
  if (show) {
    c.fillStyle = "#ffffff";
    c.textAlign = "center";
    c.textBaseline = "middle";
    c.fillText(label, left + width / 2, y + height / 2);
    c.textAlign = "left";
    c.textBaseline = "alphabetic";
  }
  return { left, right: left + width };
}
export function drawIGVGap(c, x, end, y, height, length, op) {
  c.strokeStyle = op === "N" ? IGV_SKIP : IGV_DELETION;
  c.lineWidth = op === "D" && height > 5 ? 2 : 1;
  c.setLineDash([]);
  c.beginPath();
  c.moveTo(x, y + height / 2);
  c.lineTo(end, y + height / 2);
  c.stroke();
  if (op !== "D" || length <= 1) return;
  c.font = "12px Arial,sans-serif";
  const label = String(length),
    width = c.measureText(label).width + 4;
  if (height < 8 || width >= end - x - 2) return;
  const left = (x + end - width) / 2;
  c.fillStyle = "#ffffff";
  c.fillRect(left, y, width, height);
  c.fillStyle = IGV_INSERTION;
  c.textAlign = "center";
  c.textBaseline = "middle";
  c.fillText(label, (x + end) / 2, y + height / 2);
  c.textAlign = "left";
  c.textBaseline = "alphabetic";
}
