# Phase local clone viewer — manual haplotype workspace

Open **[index.html](index.html)** directly in desktop Chrome or Edge. The viewer is a self-contained HTML file with no server, installation, uploads or runtime downloads.

Choose the variant that fits your workflow:

| Variant | Open | Approach |
| --- | --- | --- |
| Manual | [index.html](index.html) | IGV-style pileup with manual read-pair merging |
| Automatic | [automatic/index.html](automatic/index.html) | Automatic grouping of compatible read evidence |
| Assembly | [assembly/index.html](assembly/index.html) | Local assembly graph and reconstructed sequence patterns |
| Specific | [specific/index.html](specific/index.html) | Phase a supplied list of variants |

Automatic and Assembly include their original source snapshots and documentation alongside the maintained viewers.

## Load a specific locus

1. Select a coordinate-sorted short-read **BAM**, its matching **BAI or CSI**, and the alignment reference as **2bit**, or uncompressed FASTA with FAI.
2. Open the files. Enter an interval such as `chr17:7668421-7675244` using **1-based inclusive coordinates**, or find a gene in the embedded hg19/hg38 catalog and choose a transcript or individual exon. The catalog must match the BAM assembly.
3. Load the region. Reads are queried through the index for that locus; named same-contig mates are retrieved through bounded indexed lookups. There is no BAM-wide assembly or reconstruction.

The pileup uses IGV-style packing: nonoverlapping read pairs share horizontal lanes, and overlapping templates occupy separate lanes. The two ends of a pair are drawn on the same horizontal line with a connector. Packing never merges their identities or evidence. Soft clips extend from their original anchors as unaligned sequence, colored per nucleotide using desktop IGV’s palette and base-quality shading. SNPs use A/C/G/T colors; insertions use purple I-shaped markers/capped size labels, deletions use black gaps with size labels when they fit, and skipped regions use blue-gray lines. **Hover details are off by default. Right-click a read and choose Read details** for MAPQ, CIGAR, strand, flags, mate position and base information. The same menu offers **Enable hover details** / **Disable hover details**. Stored BAM sequence is already aligned in its reference orientation.

Drag horizontally to pan, use arrow/zoom buttons, or Ctrl/Command+wheel to zoom. Ordinary wheel scrolling moves vertically through the packed lanes. Large (18 px) read bars are the default; Medium and Extra large sizes are available. Panning stays within the loaded reference interval; load a different interval to change the workspace. An individual mate outside the loaded target remains available in its pair and metadata, although reference comparison is unavailable beyond the target. The initial view focuses on the first covered part of a larger target; **Fit target** shows the complete selected interval.

## Maximize the pileup

**Files & region** and **Haplotypes / VAF** independently toggle the side panels. **Compact top / Expand top** compresses the header and workspace summary while preserving selection and Undo/Redo controls. After loading a locus, both side panels are collapsed and the top is compact; the pileup fills the remaining window. Hierarchy and VAF remain available through the right-panel toggle.

## Construct haplotypes yourself

- **Click** an individual read or its pair connector to toggle that pair. Other pairs sharing its lane remain independent. **Shift-click** another pair to select a range in visual order (top to bottom, then left to right), including across scrolling. **Shift-drag** selects pairs intersecting a rectangular area. Consensus-track selection selects that group’s displayed members.
- Sort unassigned pairs by position, name, or the allele at a site selected by clicking the reference track. Sorting never changes assignments. Position mode packs densely, Name breaks same-position ties, and selected-site allele mode packs each allele band separately. **Select visible pairs** is available for batches.
- Enter an optional name and **Merge selected**. Original pairs are retained under one consensus row. Expand via the triangle or inspector to see every contributing pair. Selecting the whole haplotype plus additional compatible pairs extends it while keeping its name and identity; selecting a subset can create a separate group.
- A single contradictory callable base, deletion, or insertion **blocks the merge**. There is no majority override. Missing coverage stays unknown. Soft clips remain unaligned evidence and do not vote as reference-aligned variants.
- Select pairs and **Hide selected** to remove unhelpful rows. **Restore hidden**, **Show hidden rows**, and Undo keep hiding reversible. By default hidden pairs are excluded from consensus, coverage and VAF. Uncheck **Exclude hidden pairs from consensus / VAF** for display-only hiding.
- **Ungroup selected** detaches selected pairs; **Dissolve haplotype** returns every member to the unassigned packed pileup. Naming, merging, hiding, restoring, ungrouping and imports support Undo/Redo. Ctrl/Command+Z works outside text fields; Shift+M merges and Shift+H hides.

Nothing is automatically merged, clustered into assignments, or downsampled. You may construct any number of local haplotypes. Unassigned pairs stay independent. Changing consensus quality thresholds retains your assignments; any newly exposed conflicts are marked and consensus FASTA export stops until they are resolved.

## Hierarchy and variant fractions

The compact hierarchy compares the groups **you construct**, using average linkage on differing jointly observed SNP calls. Missing observations never become reference matches. Groups without comparable SNP evidence remain in separate trees. The similarity-cut slider changes descriptive similarity clusters without changing your groups or read assignments. Click a leaf or group label to inspect its counts and variants. Hierarchy computation is bounded to 80 constructed groups; above that, assignments/read access remain intact and the limit is explicit.

The inspector distinguishes:

- **Group pair share:** included assigned fragments divided by all included fragments in the loaded target. The denominator includes unassigned pairs and pairs from other covered blocks; it is a descriptive workspace fraction.
- **Site VAF:** fragments carrying the displayed alternate allele divided by all included fragments callable at that position, across assigned and unassigned pairs. Overlapping mates count once; uncovered, filtered and internally contradictory calls do not become reference votes.
- **Within-group ALT fraction:** alternate-carrying fragments divided by callable fragments within the selected group.

Every fraction exposes numerator and denominator. A variant shared by several groups can have a higher VAF than any one group's share. Neither group share nor site VAF measures cancer-cell fraction. Similarity does not establish tumor ancestry. Pairs with the BAM duplicate flag (0x400) are excluded from the pileup, consensus, coverage, VAF and group-share denominators by default. **Show duplicates** restores them explicitly. Unflagged pairs with identical sequence are not guessed to be PCR duplicates. All original records remain in memory. UMI family consensus is not implemented.

## Minimum variant VAF

The manual, automatic, assembly, and specific-variants viewers default to **Minimum variant VAF (%) = 2**. Alleles below 2% are hidden; exactly 2% remains included. Set the value to **0** to disable suppression. Open **Files & region** / **Evidence filters** / **Evidence settings** and change the cutoff. All variants update their rendered evidence automatically using the loaded reads; no Apply or Analyze click is needed. Read data and manual group memberships are retained.


The [specific-variants viewer](specific/index.html) has the same control under **Evidence settings**. Changing the cutoff automatically refreshes the existing results; suppressed requested variants and helpers are excluded from phase links, and the completion message reports the number suppressed.
Suppression is evaluated per allele using quality-eligible unique physical fragments at the local site, before suppression. Overlapping mates count once, internally contradictory observations do not vote, and uncovered sites are unknown. Indel absence requires a spanned insertion boundary or both deletion flanks; equivalent left-normalized events share support. Without a reference, the modal observed base supplies the comparison baseline. Structural-link annotations use a local fragment evidence fraction, an approximate measure for unresolved breakpoints rather than a resolved SV VAF.

Suppressed observations become unknown for grouping, consensus, hierarchy and cis/trans matching. They cannot cause a clash and are never rewritten as reference. The assembly viewer also excludes masked SNP/indel sequence from its k-mer paths; it does not create replacement reference sequence. Requests for suppressed alleles report the threshold and omit their highlights. Reads and original sequences remain available. Lowering the threshold recomputes the evidence; it can reveal conflicts in existing manually constructed groups without changing their memberships. Saved manual assignments include the threshold; older exports without it default to 2%.

## Save and resume

**Export assignments** writes JSON containing group names, fragment identities, hidden-pair membership and evidence settings. Restore it after loading the same BAM, reference and interval. A SHA-256 digest of the loaded reference/read evidence prevents accidental application to a different workspace. Invalid, missing or multiply assigned members are rejected atomically. Export files contain read identities and should be treated as case data.

**Export consensus FASTA** writes only observed, callable sequence blocks, with insertions incorporated and deletions removed. Gaps between covered blocks produce separate FASTA records; they are not filled from the reference or silently joined. FASTA is a partial local sequence hypothesis, not a full chromosome or cellular-clone genome. Conflicted groups cannot be exported as an arbitrary consensus.

Loading a new interval begins a fresh workspace. Export assignments to retain manual work before switching targets. **Clear session** terminates the worker and clears case data, assignments, tooltips, selected coordinates and undo history. No browser storage is used.

## Synthetic examples and validation

The built-in simple example has 52 physical fragments / 104 alignments, SNPs, insertion/deletion alleles and legal soft clips on paired reads. Request read details for the first founder pair for MAPQ 42 and 58. The second example has ten nested/sibling somatic lineages of unequal abundance plus five shared founder-only fragments. Both begin **ungrouped**, so construction remains a human decision.

`tests/fixtures/` contains generated binary BAM/BAI/CSI, 2bit and FASTA/FAI inputs. Load `synthetic:1001-1800` for the existing `synthetic`, `paired`, `record-compatibility`, `structural` and `somatic-ten` fixtures. No patient data is included or required.

```sh
npm ci
npm run build:automatic    # rebuild automatic/index.html
npm run build:assembly     # rebuild assembly/index.html
npm run build              # rebuild index.html
npm test
npm run test:browser
npm run test:browser:v2
npm run test:browser:v1
```

Browser checks use installed desktop Chrome through `playwright-core`; set `CHROME_PATH` if needed. Checks exercise real Shift-click gestures, merging/hiding, conflicts, Undo/Redo, MAPQ/soft-clip hover, read/consensus SNP geometry, shared-lane hit testing, packed layout, hierarchy, VAF denominators, duplicate filtering, right-click details and optional hover, collapsible panels, viewport sizing, scrolling/panning/zoom, assignments, binary inputs, missing SEQ, exon lookup, reset and offline CSP behavior. Real BAM testing is performed by the user.

## Implementation and bounds

New source is in `v3/`; `scripts/build.mjs` generates the top-level release and CSP hashes. Existing JavaScript BAM/reference/CIGAR readers and SNP-distance/average-linkage helpers are reused. `v2/` retains the assembly source and `src/` retains shared modules plus the original viewer source. The historical assembly algorithm review is preserved in `assembly/LITERATURE_REVIEW.md`.

Limits fail explicitly: 200 kb selected interval, 100,000 fetched records, two million stored read bases, bounded local I/O, at most 512 indexed mate windows, and a 90-second loading watchdog. There is no silent read downsampling. Read rendering is virtualized; every retained pair is reachable through scrolling. Off-contig or unmapped mates are metadata rather than fabricated aligned sequence.

Consensus defaults to MAPQ ≥20 and base quality ≥20. Secondary, QC-failed, missing/invalid-sequence and MAPQ-255 alignments do not vote. Missing/255 base quality is accepted as unknown. Missing SEQ is drawn without fake bases; inconsistent CIGAR/SEQ geometry is quarantined per record. Contig names/lengths are checked against the reference, but matching lengths do not establish sequence identity.

All dependencies and hg19/hg38 transcript annotations are embedded in the release. Runtime network requests are blocked by CSP. There are no WASM/runtime downloads. Third-party notices/build checksums are generated under `dist/`; annotation provenance is in `references/manifest.json`. Case inputs are ignored by `.gitignore`, with intentional exceptions for generated synthetic fixtures.

The packing/display behavior follows the [IGV alignment guide](https://igv.org/doc/desktop/UserGuide/tracks/alignments/viewing_alignments_basics/) and [IGV paired-end view](https://igv.org/doc/desktop/UserGuide/tracks/alignments/paired_end_alignments/). Phase uses its own local Canvas renderer and packing implementation; no IGV server or runtime download is required.

## Desktop IGV base and indel styling

The renderer follows desktop IGV light-theme defaults at commit `811bf7042e0b05d124ff421ee7100eec4789d179`: [palette/default preferences](https://github.com/igvteam/igv/blob/811bf7042e0b05d124ff421ee7100eec4789d179/src/main/resources/preferences.tab), [alignment/indel rendering](https://github.com/igvteam/igv/blob/811bf7042e0b05d124ff421ee7100eec4789d179/src/main/java/org/igv/alignment/AlignmentRenderer.java), and [base rendering/quality shading](https://github.com/igvteam/igv/blob/811bf7042e0b05d124ff421ee7100eec4789d179/src/main/java/org/igv/alignment/BaseRenderer.java).

A = `#009600`, C = `#0000ff`, G = `#d17105`, T = `#ff0000`, N = `#808080`. Soft clips use those colors per base, not a uniform purple overlay. Bases are colored blocks at overview zoom and colored letters on the gray read background once base width permits. Quality shading uses IGV's default Q5/Q20 bounds and tenths quantization. Insertions use `#7618dc` capped I markers and centered size labels when they fit; deletions are black, with white-backed size labels; skipped regions are `#96b8c8`. Manual haplotype identification colors remain distinct from nucleotide colors. These match standard desktop alignment styling, not every optional IGV mode (e.g. base modifications or user-customized themes).

Reference bases and read bases share one horizontal coordinate axis and the same centered letter renderer, including when a vertical scrollbar reduces alignment-track width. Insertion markers on original reads and merged haplotypes show inserted sequence on hover or click, even with general read-hover details disabled. A plain insertion-marker click inspects sequence without changing read selection; Shift-click retains range-selection behavior.
