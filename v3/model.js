import { applyVafFilter, hiddenBase, hiddenEvent } from "../src/vaf-filter.js";
import {
  prepareRecord,
  alignmentRuns,
  hasSequence,
} from "../src/alignment-geometry.js";

export const fragmentKey = (r) =>
  JSON.stringify([r.sample || "", r.rg || "", r.name]);
const add = (map, pos, value) => {
  if (!map.has(pos)) map.set(pos, new Set());
  map.get(pos).add(value);
};
const sameMembers = (a, b) => a.length === b.size && a.every((id) => b.has(id));
export function recordEvidence(r, region, { baseq = 20, mapq = 20 } = {}) {
  const bases = new Map(),
    insertions = new Map();
  if (
    !hasSequence(r) ||
    !r.alignmentUsable ||
    r.flags & (4 | 256 | 512) ||
    r.mq === 255 ||
    r.mq < mapq
  )
    return { bases, insertions };
  const good = (q) =>
    q >= 0 &&
    q < r.seq.length &&
    (!r.qual?.length || r.qual[q] === 255 || r.qual[q] >= baseq);
  const runs = r.runs || alignmentRuns(r);
  for (const run of runs) {
    const { op, length: n, refStart: p, queryStart: q } = run;
    if ("M=X".includes(op))
      for (let j = 0; j < n; j++) {
        const pos = p + j,
          b = r.seq[q + j];
        if (
          pos < region.start ||
          pos >= region.end ||
          !good(q + j) ||
          hiddenBase(r, q + j) ||
          !/[ACGT]/.test(b)
        )
          continue;
        bases.set(pos, b);
        // An observed contiguous reference step establishes absence of an
        // insertion. Coverage on one side alone never implies a reference call.
        if (j > 0 && good(q + j - 1)) insertions.set(pos, "");
      }
    else if (
      op === "I" &&
      !hiddenEvent(r, op, p) &&
      p >= region.start &&
      p < region.end &&
      Array.from({ length: n }, (_, j) => good(q + j)).every(Boolean)
    )
      insertions.set(p, r.seq.slice(q, q + n));
    else if (op === "D" && !hiddenEvent(r, op, p) && good(q - 1) && good(q))
      for (
        let pos = Math.max(p, region.start);
        pos < Math.min(p + n, region.end);
        pos++
      )
        bases.set(pos, "-");
  }
  return { bases, insertions };
}
export class ManualWorkspace {
  constructor(data) {
    this.region = data.region;
    this.reference = data.reference;
    this.dataset = data.dataset;
    this.reads = data.reads;
    this.settings = {
      mapq: 20,
      baseq: 20,
      excludeHidden: true,
      includeDuplicates: false,
      minVaf: 0.02,
    };
    this.hidden = new Set();
    this.groups = [];
    this.selection = new Set();
    this.expanded = new Set();
    this.nextID = 1;
    this.undoStack = [];
    this.redoStack = [];
    this.cache = new Map();
    this.groupCache = new Map();
    const byKey = new Map();
    this.fragments = [];
    this.readById = new Map();
    for (const r of this.reads) {
      prepareRecord(r);
      r.runs = alignmentRuns(r);
      this.readById.set(r.id, r);
      const key = fragmentKey(r);
      if (!byKey.has(key)) {
        const f = {
          id: this.fragments.length,
          key,
          name: r.name,
          readIds: [],
          duplicate: false,
          start: Infinity,
          end: 0,
        };
        byKey.set(key, f);
        this.fragments.push(f);
      }
      const f = byKey.get(key);
      f.readIds.push(r.id);
      if (r.flags & 1024) f.duplicate = true;
      if (r.alignmentUsable) {
        f.start = Math.min(f.start, r.start);
        f.end = Math.max(f.end, r.end);
      }
    }
    this.fragmentByKey = byKey;
    this.owner = new Map();
    this.updateVaf();
  }
  updateVaf() {
    const active = new Set(
      this.activeIDs(this.fragments.map((f) => f.id)).flatMap(
        (id) => this.fragments[id].readIds,
      ),
    );
    this.vafFilter = applyVafFilter(this.reads, this.region, this.reference, {
      ...this.settings,
      eligible: (r) => active.has(r.id),
    });
  }
  setSettings(settings) {
    if (
      settings.minVaf !== undefined &&
      (!Number.isFinite(settings.minVaf) ||
        settings.minVaf < 0 ||
        settings.minVaf > 1)
    )
      throw Error("Invalid minimum VAF");
    if (
      settings.includeDuplicates !== undefined &&
      typeof settings.includeDuplicates !== "boolean"
    )
      throw Error("Invalid duplicate filter");
    for (const k of ["mapq", "baseq"])
      if (
        !Number.isInteger(settings[k]) ||
        settings[k] < 0 ||
        settings[k] > (k === "mapq" ? 254 : 93)
      )
        throw Error("Invalid consensus quality thresholds");
    this.settings = { ...this.settings, ...settings };
    this.cache.clear();
    this.groupCache.clear();
    this.updateVaf();
  }
  snapshot() {
    return {
      groups: this.groups.map((g) => ({ ...g, members: [...g.members] })),
      hidden: [...this.hidden],
      settings: { ...this.settings },
      nextID: this.nextID,
    };
  }
  applySnapshot(s) {
    this.groups = s.groups.map((g) => ({ ...g, members: [...g.members] }));
    this.hidden = new Set(s.hidden || []);
    this.setSettings(s.settings || this.settings);
    this.nextID = s.nextID;
    this.owner.clear();
    for (const g of this.groups)
      for (const id of g.members) this.owner.set(id, g.id);
    this.selection.clear();
    this.groupCache.clear();
    for (const id of this.expanded)
      if (!this.groups.some((g) => g.id === id)) this.expanded.delete(id);
  }
  commit(change) {
    const before = this.snapshot();
    change();
    this.undoStack.push(before);
    if (this.undoStack.length > 40) this.undoStack.shift();
    this.redoStack = [];
    this.applySnapshot(this.snapshot());
  }
  undo() {
    if (!this.undoStack.length) return;
    const s = this.undoStack.pop();
    this.redoStack.push(this.snapshot());
    this.applySnapshot(s);
  }
  redo() {
    if (!this.redoStack.length) return;
    const s = this.redoStack.pop();
    this.undoStack.push(this.snapshot());
    this.applySnapshot(s);
  }
  evidence(id) {
    if (this.cache.has(id)) return this.cache.get(id);
    const f = this.fragments[id],
      bases = new Map(),
      insertions = new Map();
    for (const rid of f.readIds) {
      const e = recordEvidence(
        this.readById.get(rid),
        this.region,
        this.settings,
      );
      for (const [p, b] of e.bases) add(bases, p, b);
      for (const [p, s] of e.insertions) add(insertions, p, s);
    }
    const e = { bases, insertions };
    this.cache.set(id, e);
    return e;
  }
  activeIDs(ids) {
    return ids.filter(
      (id) =>
        (this.settings.includeDuplicates || !this.fragments[id].duplicate) &&
        (!this.settings.excludeHidden || !this.hidden.has(id)),
    );
  }
  hideSelection() {
    if (!this.selection.size) throw Error("Select pairs to hide");
    const ids = [...this.selection];
    this.commit(() => {
      for (const id of ids) this.hidden.add(id);
    });
  }
  restoreSelection() {
    const ids = [...this.selection];
    this.commit(() => {
      for (const id of ids) this.hidden.delete(id);
    });
  }
  restoreAll() {
    if (!this.hidden.size) return;
    this.commit(() => this.hidden.clear());
  }
  summarize(ids) {
    ids = this.activeIDs(ids);
    const bases = new Map(),
      insertions = new Map(),
      clips = new Map();
    let unknown = 0;
    const accumulate = (out, p, value, id) => {
      if (!out.has(p)) out.set(p, new Map());
      const byValue = out.get(p);
      if (!byValue.has(value)) byValue.set(value, new Set());
      byValue.get(value).add(id);
    };
    for (const id of ids) {
      const f = this.fragments[id],
        e = this.evidence(id);
      if (!e.bases.size && !e.insertions.size) unknown++;
      for (const [p, bs] of e.bases)
        for (const b of bs) accumulate(bases, p, b, id);
      for (const [p, ss] of e.insertions)
        for (const s of ss) accumulate(insertions, p, s, id);
      for (const rid of f.readIds) {
        const r = this.readById.get(rid);
        if (!hasSequence(r)) continue;
        for (const run of r.runs)
          if (run.op === "S") {
            const side = run.queryStart === 0 ? "left" : "right",
              seq = r.seq.slice(run.queryStart, run.queryStart + run.length),
              key = JSON.stringify([run.refStart, side, seq]);
            if (!clips.has(key))
              clips.set(key, {
                pos: run.refStart,
                side,
                seq,
                fragments: new Set(),
              });
            clips.get(key).fragments.add(id);
          }
      }
    }
    const calls = (m) =>
      [...m]
        .sort((a, b) => a[0] - b[0])
        .map(([pos, values]) => ({
          pos,
          values: [...values].map(([value, fs]) => ({
            value,
            support: fs.size,
          })),
          conflict: values.size > 1,
          total: new Set([...values.values()].flatMap((s) => [...s])).size,
        }));
    const bs = calls(bases),
      ins = calls(insertions).filter(
        (c) => c.conflict || c.values[0].value.length,
      ),
      conflicts = [
        ...bs
          .filter((c) => c.conflict)
          .map((c) => ({ ...c, kind: "base/deletion" })),
        ...ins
          .filter((c) => c.conflict)
          .map((c) => ({ ...c, kind: "insertion" })),
      ];
    return {
      bases: bs,
      insertions: ins,
      conflicts,
      unknown,
      clips: [...clips.values()].map((c) => ({
        ...c,
        support: c.fragments.size,
        fragments: undefined,
      })),
      fragments: ids.length,
    };
  }
  groupSummary(g) {
    const key = g.id + ":" + g.members.join(",");
    if (!this.groupCache.has(key))
      this.groupCache.set(key, this.summarize(g.members));
    return this.groupCache.get(key);
  }
  merge(name = "") {
    const ids = [...this.selection].sort((a, b) => a - b);
    if (!ids.length) throw Error("Select at least one read pair first");
    const summary = this.summarize(ids);
    if (summary.conflicts.length) {
      const c = summary.conflicts[0];
      throw Error(
        `Cannot merge: contradictory ${c.kind} at ${this.region.chr}:${c.pos + 1} (${c.values.map((v) => v.value || "no insertion").join(" / ")}). ${summary.conflicts.length} conflicting site(s).`,
      );
    }
    const existing = this.groups.find((g) =>
      sameMembers(g.members, this.selection),
    );
    if (existing) {
      if (name.trim()) this.rename(existing.id, name);
      this.selection.clear();
      return existing.id;
    }
    // Iterative manual extension keeps the existing group's identity, color
    // and name. Selecting a subset still creates a separate new hypothesis.
    const touched = this.groups.filter((g) =>
      g.members.some((id) => this.selection.has(id)),
    );
    if (
      touched.length === 1 &&
      touched[0].members.every((id) => this.selection.has(id))
    ) {
      const old = touched[0],
        label = (name.trim() || old.name).slice(0, 80);
      this.commit(() => {
        this.groups = this.groups.map((g) =>
          g.id === old.id ? { ...g, name: label, members: ids } : g,
        );
      });
      return old.id;
    }
    const id = this.nextID,
      label = (name.trim() || `H${id}`).slice(0, 80),
      selected = new Set(ids);
    this.commit(() => {
      this.groups = this.groups
        .map((g) => ({
          ...g,
          members: g.members.filter((id) => !selected.has(id)),
        }))
        .filter((g) => g.members.length);
      this.groups.push({ id, name: label, members: ids });
      this.nextID++;
    });
    return id;
  }
  detach() {
    if (!this.selection.size) throw Error("Select pairs to ungroup");
    const ids = new Set(this.selection);
    this.commit(() => {
      this.groups = this.groups
        .map((g) => ({ ...g, members: g.members.filter((id) => !ids.has(id)) }))
        .filter((g) => g.members.length);
    });
  }
  dissolve(id) {
    if (!this.groups.some((g) => g.id === id)) return;
    this.commit(() => {
      this.groups = this.groups.filter((g) => g.id !== id);
    });
  }
  rename(id, name) {
    name = name.trim().slice(0, 80);
    if (!name) throw Error("Enter a haplotype name");
    if (!this.groups.some((g) => g.id === id))
      throw Error("Haplotype no longer exists");
    this.commit(() => {
      this.groups = this.groups.map((g) => (g.id === id ? { ...g, name } : g));
    });
  }
  rows(order = "position", site = null, showHidden = false) {
    const rows = [];
    for (const g of this.groups) {
      rows.push({
        key: "g:" + g.id,
        kind: "group",
        group: g,
        members: g.members.filter(
          (id) =>
            (this.settings.includeDuplicates ||
              !this.fragments[id].duplicate) &&
            (showHidden || !this.hidden.has(id)),
        ),
      });
      if (this.expanded.has(g.id))
        for (const id of g.members)
          if (
            (this.settings.includeDuplicates ||
              !this.fragments[id].duplicate) &&
            (showHidden || !this.hidden.has(id))
          )
            rows.push({
              key: "c:" + id,
              kind: "child",
              fragment: this.fragments[id],
              members: [id],
              parent: g.id,
            });
    }
    const free = this.fragments.filter(
      (f) =>
        !this.owner.has(f.id) &&
        (this.settings.includeDuplicates || !f.duplicate) &&
        (showHidden || !this.hidden.has(f.id)),
    );
    const allele = (id) => {
      if (site === null) return "";
      const v = this.evidence(id).bases.get(site);
      return v?.size === 1 ? [...v][0] : v?.size > 1 ? "!" : "?";
    };
    free.sort(
      order === "name"
        ? (a, b) => a.name.localeCompare(b.name)
        : order === "allele"
          ? (a, b) =>
              allele(a.id).localeCompare(allele(b.id)) || a.start - b.start
          : (a, b) => a.start - b.start || a.name.localeCompare(b.name),
    );
    for (const f of free)
      rows.push({
        key: "f:" + f.id,
        kind: "fragment",
        fragment: f,
        members: [f.id],
      });
    return rows;
  }
  coverage() {
    const n = this.region.end - this.region.start,
      delta = new Int32Array(n + 1);
    for (const f of this.fragments) {
      if (!this.settings.includeDuplicates && f.duplicate) continue;
      if (this.settings.excludeHidden && this.hidden.has(f.id)) continue;
      const intervals = [];
      for (const id of f.readIds) {
        const r = this.readById.get(id);
        if (r.flags & (4 | 256 | 2048)) continue;
        for (const run of r.runs)
          if ("M=X".includes(run.op)) {
            const a = Math.max(this.region.start, run.refStart),
              b = Math.min(this.region.end, run.refStart + run.length);
            if (b > a) intervals.push([a, b]);
          }
      }
      intervals.sort((a, b) => a[0] - b[0]);
      let merged = [];
      for (const [a, b] of intervals) {
        const last = merged.at(-1);
        if (last && a <= last[1]) last[1] = Math.max(last[1], b);
        else merged.push([a, b]);
      }
      for (const [a, b] of merged) {
        delta[a - this.region.start]++;
        delta[b - this.region.start]--;
      }
    }
    let depth = 0;
    for (let i = 0; i < n; i++) {
      depth += delta[i];
      delta[i] = depth;
    }
    return delta.subarray(0, n);
  }
  exportAssignments() {
    return JSON.stringify(
      {
        format: "phase-manual-assignments",
        version: 1,
        dataset: this.dataset,
        region: this.region,
        settings: this.settings,
        hidden: [...this.hidden].map((id) => this.fragments[id].key),
        groups: this.groups.map((g) => ({
          name: g.name,
          members: g.members.map((id) => this.fragments[id].key),
        })),
      },
      null,
      2,
    );
  }
  importAssignments(text) {
    if (text.length > 16000000) throw Error("Assignment file exceeds 16 MB");
    const d = JSON.parse(text);
    if (
      d.format !== "phase-manual-assignments" ||
      d.version !== 1 ||
      d.dataset !== this.dataset ||
      JSON.stringify(d.region) !== JSON.stringify(this.region)
    )
      throw Error(
        "Assignments do not match this BAM/reference interval. Load the same input and coordinates first.",
      );
    if (!Array.isArray(d.groups) || d.groups.length > this.fragments.length)
      throw Error("Invalid assignment groups");
    const seen = new Set(),
      groups = d.groups.map((g, i) => {
        if (
          typeof g.name !== "string" ||
          !g.name.trim() ||
          g.name.length > 80 ||
          !Array.isArray(g.members) ||
          !g.members.length
        )
          throw Error("Invalid haplotype entry");
        const members = g.members.map((key) => {
          const f = this.fragmentByKey.get(key);
          if (!f || seen.has(f.id))
            throw Error("Missing or multiply assigned read pair");
          seen.add(f.id);
          return f.id;
        });
        return { id: i + 1, name: g.name, members };
      });
    const s = d.settings;
    if (
      !s ||
      (s.minVaf !== undefined &&
        (!Number.isFinite(s.minVaf) || s.minVaf < 0 || s.minVaf > 1)) ||
      !Number.isInteger(s.mapq) ||
      s.mapq < 0 ||
      s.mapq > 254 ||
      !Number.isInteger(s.baseq) ||
      s.baseq < 0 ||
      s.baseq > 93 ||
      typeof s.excludeHidden !== "boolean" ||
      (s.includeDuplicates !== undefined &&
        typeof s.includeDuplicates !== "boolean")
    )
      throw Error("Invalid consensus thresholds in assignment file");
    if (!Array.isArray(d.hidden) || d.hidden.length > this.fragments.length)
      throw Error("Invalid hidden pair list");
    const hidden = new Set(
      d.hidden.map((key) => {
        const f = this.fragmentByKey.get(key);
        if (!f) throw Error("Hidden pair not in this interval");
        return f.id;
      }),
    );
    this.commit(() => {
      this.groups = groups;
      this.hidden = hidden;
      this.nextID = groups.length + 1;
      this.setSettings({
        ...s,
        includeDuplicates: s.includeDuplicates ?? false,
        minVaf: s.minVaf ?? 0.02,
      });
    });
  }
  exportFasta() {
    const out = [];
    for (const g of this.groups) {
      const s = this.groupSummary(g);
      if (s.conflicts.length)
        throw Error(
          `Resolve conflicting evidence in ${g.name} before exporting consensus`,
        );
      const byPos = new Map(s.bases.map((c) => [c.pos, c.values[0].value])),
        ins = new Map(s.insertions.map((c) => [c.pos, c.values[0].value]));
      const ps = [...byPos.keys()].sort((a, b) => a - b);
      let start = null,
        last = null,
        seq = "";
      const flush = () => {
        if (start !== null && seq)
          out.push(
            `>${g.name.replace(/[^A-Za-z0-9_.-]/g, "_")}|${this.region.chr}:${start + 1}-${last + 1}|pairs=${g.members.length}|observed_block\n${seq.match(/.{1,80}/g).join("\n")}`,
          );
      };
      for (const p of ps) {
        if (last === null || p !== last + 1) {
          flush();
          start = p;
          seq = "";
        } else if (ins.has(p)) seq += ins.get(p);
        const b = byPos.get(p);
        if (b !== "-") seq += b;
        last = p;
      }
      flush();
    }
    if (!out.length)
      throw Error("Create a haplotype with callable sequence first");
    return out.join("\n") + "\n";
  }
}
