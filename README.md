# Phase local clone viewer — manual haplotype workspace

Open **[index.html](index.html)** directly in desktop Chrome or Edge. The viewer is a self-contained HTML file with no server, installation, uploads or runtime downloads.

The assembly-first version is preserved at **[archive2/index.html](archive2/index.html)**, with its source snapshot, documentation, literature notes and checksums in the same folder. The original viewer remains at **[archive/index.html](archive/index.html)**.

## Load a specific locus

1. Select a coordinate-sorted short-read **BAM**, its matching **BAI or CSI**, and the alignment reference as **2bit**, or uncompressed FASTA with FAI.
2. Open the files. Enter an interval such as `chr17:7668421-7675244` using **1-based inclusive coordinates**, or find a gene in the embedded hg19/hg38 catalog and choose a transcript or individual exon. The catalog must match the BAM assembly.
3. Load the region. Reads are queried through the index for that locus; named same-contig mates are retrieved through bounded indexed lookups. There is no BAM-wide assembly or reconstruction.

Each physical read pair starts as its own row. The two ends share a coordinate axis and are connected across the insert. Soft clips extend from their original anchors as lighter, unaligned sequence. SNPs use A/C/G/T colors; insertions are purple triangles; deletions and skipped regions are drawn separately. Hover reads for **MAPQ**, CIGAR, strand, flags, mate position and base information. Stored BAM sequence is already aligned in its reference orientation.

Drag horizontally to pan, use arrow/zoom buttons, or Ctrl/Command+wheel to zoom. Ordinary wheel scrolling moves vertically through all rows. Panning stays within the loaded reference interval; load a different interval to change the workspace. An individual mate outside the loaded target remains available in its pair and metadata, although reference comparison is unavailable beyond the target. The initial view focuses on the first covered part of a larger target; **Fit target** shows the complete selected interval.

## Construct haplotypes yourself

- **Click** a row to toggle its whole pair. **Shift-click** another row to select the inclusive range in the current order, including rows reached by scrolling. Group-row selection selects its displayed members; expanded child rows select individual pairs.
- Sort unassigned pairs by position, name, or the allele at a site selected by clicking the reference track. Sorting never changes assignments. **Select visible rows** is available for batches.
- Enter an optional name and **Merge selected**. Original pairs are retained under one consensus row. Expand via the triangle or inspector to see every contributing pair. Selecting the whole haplotype plus additional compatible pairs extends it while keeping its name and identity; selecting a subset can create a separate group.
- A single contradictory callable base, deletion, or insertion **blocks the merge**. There is no majority override. Missing coverage stays unknown. Soft clips remain unaligned evidence and do not vote as reference-aligned variants.
- Select pairs and **Hide selected** to remove unhelpful rows. **Restore hidden**, **Show hidden rows**, and Undo keep hiding reversible. By default hidden pairs are excluded from consensus, coverage and VAF. Uncheck **Exclude hidden pairs from consensus / VAF** for display-only hiding.
- **Ungroup selected** detaches selected pairs; **Dissolve haplotype** returns every member to individual rows. Naming, merging, hiding, restoring, ungrouping and imports support Undo/Redo. Ctrl/Command+Z works outside text fields; Shift+M merges and Shift+H hides.

Nothing is automatically merged, clustered into assignments, or downsampled. You may construct any number of local haplotypes. Unassigned pairs stay independent. Changing consensus quality thresholds retains your assignments; any newly exposed conflicts are marked and consensus FASTA export stops until they are resolved.

## Hierarchy and variant fractions

The compact hierarchy compares the groups **you construct**, using average linkage on differing jointly observed SNP calls. Missing observations never become reference matches. Groups without comparable SNP evidence remain in separate trees. The similarity-cut slider changes descriptive similarity clusters without changing your groups or read assignments. Click a leaf or group label to inspect its counts and variants. Hierarchy computation is bounded to 80 constructed groups; above that, assignments/read access remain intact and the limit is explicit.

The inspector distinguishes:

- **Group pair share:** included assigned fragments divided by all included fragments in the loaded target. The denominator includes unassigned pairs and pairs from other covered blocks; it is a descriptive workspace fraction.
- **Site VAF:** fragments carrying the displayed alternate allele divided by all included fragments callable at that position, across assigned and unassigned pairs. Overlapping mates count once; uncovered, filtered and internally contradictory calls do not become reference votes.
- **Within-group ALT fraction:** alternate-carrying fragments divided by callable fragments within the selected group.

Every fraction exposes numerator and denominator. A variant shared by several groups can have a higher VAF than any one group's share. Neither group share nor site VAF measures cancer-cell fraction. Similarity does not establish tumor ancestry. Duplicate-flagged pairs remain available and, when manually selected, contribute as distinct query-name fragments; UMI family consensus is not implemented.

## Save and resume

**Export assignments** writes JSON containing group names, fragment identities, hidden-pair membership and evidence settings. Restore it after loading the same BAM, reference and interval. A SHA-256 digest of the loaded reference/read evidence prevents accidental application to a different workspace. Invalid, missing or multiply assigned members are rejected atomically. Export files contain read identities and should be treated as case data.

**Export consensus FASTA** writes only observed, callable sequence blocks, with insertions incorporated and deletions removed. Gaps between covered blocks produce separate FASTA records; they are not filled from the reference or silently joined. FASTA is a partial local sequence hypothesis, not a full chromosome or cellular-clone genome. Conflicted groups cannot be exported as an arbitrary consensus.

Loading a new interval begins a fresh workspace. Export assignments to retain manual work before switching targets. **Clear session** terminates the worker and clears case data, assignments, tooltips, selected coordinates and undo history. No browser storage is used.

## Synthetic examples and validation

The built-in simple example has 52 physical fragments / 104 alignments, SNPs, insertion/deletion alleles and legal soft clips on paired reads. Hover the first founder pair for MAPQ 42 and 58. The second example has ten nested/sibling somatic lineages of unequal abundance plus five shared founder-only fragments. Both begin **ungrouped**, so construction remains a human decision.

`tests/fixtures/` contains generated binary BAM/BAI/CSI, 2bit and FASTA/FAI inputs. Load `synthetic:1001-1800` for the existing `synthetic`, `paired`, `record-compatibility`, `structural` and `somatic-ten` fixtures. No patient data is included or required.

```sh
npm ci
npm run build
npm test
npm run test:browser
npm run test:browser:v2
npm run test:browser:v1
```

Browser checks use installed desktop Chrome through `playwright-core`; set `CHROME_PATH` if needed. Checks exercise real Shift-click gestures, merging/hiding, conflicts, Undo/Redo, MAPQ/soft-clip hover, read/consensus SNP geometry, hierarchy, VAF denominators, scrolling/panning/zoom, assignments/FASTA, binary inputs, missing SEQ, exon lookup, reset and offline CSP behavior. Real BAM testing is performed by the user.

## Implementation and bounds

New source is in `v3/`; `scripts/build.mjs` generates the top-level release and CSP hashes. Existing JavaScript BAM/reference/CIGAR readers and SNP-distance/average-linkage helpers are reused. `v2/` retains the assembly source and `src/` retains shared modules plus the original viewer source. The historical assembly algorithm review is preserved in `archive2/LITERATURE_REVIEW.md`.

Limits fail explicitly: 200 kb selected interval, 100,000 fetched records, two million stored read bases, bounded local I/O, at most 512 indexed mate windows, and a 90-second loading watchdog. There is no silent read downsampling. Read rendering is virtualized; every retained pair is reachable through scrolling. Off-contig or unmapped mates are metadata rather than fabricated aligned sequence.

Consensus defaults to MAPQ ≥20 and base quality ≥20. Secondary, QC-failed, missing/invalid-sequence and MAPQ-255 alignments do not vote. Missing/255 base quality is accepted as unknown. Missing SEQ is drawn without fake bases; inconsistent CIGAR/SEQ geometry is quarantined per record. Contig names/lengths are checked against the reference, but matching lengths do not establish sequence identity.

All dependencies and hg19/hg38 transcript annotations are embedded in the release. Runtime network requests are blocked by CSP. There are no WASM/runtime downloads. Third-party notices/build checksums are generated under `dist/`; annotation provenance is in `references/manifest.json`. Case inputs are ignored by `.gitignore`, with intentional exceptions for generated synthetic fixtures.
