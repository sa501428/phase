# Phase local clone viewer — assembly version

Open **[index.html](index.html)** directly in desktop Chrome or Edge. It is a self-contained HTML file: no server, installation, uploads or runtime downloads. The previous viewer remains directly usable at **[archive/index.html](archive/index.html)**. Its frozen source snapshot and original documentation are under `archive/v1-2026-11-10/`.

## Work on a selected locus

1. Choose a coordinate-sorted short-read BAM and matching BAI or CSI. Optionally choose the matching reference as `.2bit`, or uncompressed FASTA with `.fai`.
2. Open the BAM. Enter a contig and **1-based inclusive** start/end, or find a gene in the embedded hg19/hg38 catalog, choose a transcript and optionally one exon. An interval spanning a few exons also works; disconnected evidence stays in separate local blocks. The chosen catalog must match the BAM assembly.
3. Assemble the selected region. The reference is optional for sequence assembly; BAM coordinates still define the local target. Off-target mate sequences do not seed another assembly.
4. Explore the compact graph. Drag left/right/up/down to pan; wheel or buttons zoom. Select a sequence segment or local haplotype pattern to inspect its supporting sequence and fragments. The graph layout is topological, not a genomic coordinate axis.
5. Expand paired read evidence when needed. Both ends share one row/axis. Scroll vertically through all fragments; drag a row horizontally to pan, or use Ctrl/Command+wheel to zoom. Rendering is virtualized; records are not downsampled. **Show all input records** also exposes filtered/malformed/missing-sequence records.
6. Enter variants, one per line, and check linkage. Genomic SNVs such as `chr17:7674220G>A` work without reference with REF explicitly unverified. With reference, the existing resolver supports genomic HGVS, transcript `c.` entries, selected-transcript bare `c.` entries, substitutions, insertions, deletions, duplications and delins. Only the selected interval is evaluated.
7. Clear session when done. All working case data is held in RAM. GFA export is explicit and contains assembled case sequence.

## What the graph and counts mean

Read sequences construct a strand-symmetric de Bruijn graph, compacted into unitigs. Low-support branches and cycles remain visible. Actual read traversals define bubble allele sequences, and physical paired fragments link allele features into local evidence blocks. A single weak bridge does not join blocks at the default support threshold.

There is **no fixed ploidy or requested number of haplotypes**. Suggestions are compatible local sequence patterns, not a measured number of cancer cell clones. A contradictory callable SNP vetoes grouping, including singleton SNPs outside the bubble list. Missing/low-quality bases are unknown. Partial evidence compatible with multiple patterns remains shared; weak patterns are inspectable separately. This is a conservative heuristic, not a globally optimal or statistically calibrated mixture model.

Pattern share is exclusive assigned fragments divided by all informative fragments in the block, including shared and contradictory evidence. Bubble allele fraction uses noncontradictory fragments callable at that bubble. Both denominators are exposed; neither measures tumor-cell fraction. Duplicate flags, PCR/amplicon bias, copy number and unequal coverage affect them. UMI consensus is not implemented.

Reference annotation maps existing graph sequence to the supplied interval; it never creates a reference-supported assembly path. Matching contig names/lengths do not verify sequence identity. Original CIGAR geometry determines read markers, including soft clips, insertions and deletions. Missing SEQ is displayed without fabricated bases; malformed records are quarantined without aborting the locus.

Requested-variant cis evidence means the alternate alleles occurred on the same physical fragment. Separate allele patterns are compatible with trans, but cannot establish trans for an unrestricted somatic mixture. The inspector's closest-pattern comparison describes shared/different observed bubble alleles; it does not infer ancestry.

The structural panel reports local SA/supplementary and mate clues for possible translocations, inversion-like alignments, and other regional links. Repeat/duplication loops and deletion/insertion branches may also appear in the graph. These are evidence for review, not exhaustive structural variant calls or resolved remote breakpoints. No partner-region assembly is attempted.

## Synthetic examples and checks

**Three patterns + deletion block** contains 52 paired fragments (104 records): three linked patterns of 12 fragments each and a separate two-allele deletion block of 8 fragments per allele.

**Ten branched somatic haplotypes** contains ten sequence lineages with nested founder/private SNPs and unequal abundances of 36, 24, 18, 12, 6, 30, 21, 15, 9 and 3 paired fragments. Five founder-only partial observations stay shared among five compatible patterns. The expected result is ten supported local patterns, not a forced binary partition. The lineage tree is synthetic ground truth, never an assembly input or a claimed inferred tumor tree.

`tests/fixtures/` contains only generated synthetic binary BAM/BAI/CSI, FASTA/FAI and 2bit data. Existing `synthetic`, `paired`, `record-compatibility`, `structural`, and new `somatic-ten` fixtures use `synthetic:1001-1800`. The latter includes all ten unequal lineages plus five shared founder-only fragments in a multi-block BAM. Unlike the old viewer, a deletion present in all reads may be a single graph path rather than an alternate bubble; optional reference and read CIGARs expose the difference.

```sh
npm ci
npm run build
npm test
npm run test:browser
npm run test:browser:v1
```

Browser checks use installed desktop Chrome via `playwright-core`; set `CHROME_PATH` if needed. They open the release through `file://`, use synthetic inputs, and verify that there are no network requests or CSP errors. Clinical BAM testing remains the user's responsibility; no private BAMs are included or required.

## Bounds and implementation

Analysis is explicitly local: maximum 200 kb selected interval, 100,000 fetched records, two million stored read bases, 200,000 oriented k-mers and 400,000 oriented edges. The worker has a 90-second watchdog. Exceeding a limit fails with a request to narrow the target; it never silently downsamples. Named same-contig mates are retrieved through at most 512 indexed windows. Cross-contig/unmapped/unavailable mates remain metadata. Local I/O has its own bounded budgets.

Defaults are k=31, MAPQ ≥20, base quality ≥20, and support ≥3 distinct fragments. Secondary, QC-failed, duplicate-flagged (unless enabled), unavailable-sequence and malformed alignments do not vote. MAPQ 255 is excluded. Missing/255 base qualities are unknown and currently accepted, which is reported in settings. These exploration thresholds are not validated limits of detection. Weak graph paths are retained. Bubble searches are bounded to nearby reconvergence (512 visited nodes per branch, at most 128 resolved bubbles); complex structures remain unresolved and counted.

New source is in `v2/`; `scripts/build.mjs` generates the top-level `index.html` and its CSP hashes. Selected low-level BAM/reference/CIGAR/variant/structural modules are reused from `src/`. The earlier viewer's pileup controls, VCF/BED flows, event realignment and dendrogram remain in the frozen archive. V2 emphasizes assembly graph evidence with an inspector rather than bottom tables.

The implementation is JavaScript, including the bounded BGZF inflater. `@gmod/bam` and `pako` are bundled; there are no WASM or runtime network dependencies. Embedded hg19/hg38 transcript provenance is in `references/manifest.json`. Third-party notices and build checksums are generated under `dist/`. Case files are ignored by `.gitignore`; only intentional synthetic fixture exceptions are tracked.

See **[LITERATURE_REVIEW.md](LITERATURE_REVIEW.md)** for the relevant Lancet, SAVAGE, WhatsHap polyphase and Haploflow ideas, their implementation consequences, and the substantial parts of those algorithms that Phase does not implement or claim to reproduce.
