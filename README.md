# Phase local clone viewer

Open **[index.html](index.html)** directly in desktop Chrome or Edge. No installation, server, or internet connection is required. The release is a single self-contained HTML file at the top level of `phase/`. `npm run build` regenerates this file; `src/index.html` is its source template. Firefox has not been verified.

## Use

1. Select a coordinate-sorted BAM and its matching BAI or CSI.
2. Select the exact alignment reference as `.2bit`, or uncompressed FASTA with its `.fai`. Confirm reference identity. Contig names/lengths are checked, but a matching length does not establish sequence identity.
3. Optionally select a plain VCF or a BED file. BED coordinates are zero-based, half-open; choose capture targets for display or primer intervals for base-quality masking.
4. Open the workspace. Choose the hg19/hg38 gene catalog, enter an official gene symbol, and select a locus/transcript, or enter a one-based inclusive interval such as `chr17:7571720-7590868`.
5. Analyze. Drag the pileup horizontally to pan (zoom in first when viewing the full interval), use the pan/zoom controls, change read sorting/layers, inspect phase blocks, select variants for joint allele counts, or open a structural event. Changed evidence filters take effect on the next analysis.
6. Use **Clear session** before opening another case. Reload also starts empty.

**Explore synthetic example** runs the same evidence engine on generated reads without opening any files. To test actual binary input, use the four matching files in `tests/fixtures`: `synthetic.bam`, `synthetic.bai`, `synthetic.fa`, and `synthetic.fa.fai`. Analyze `synthetic:1001-1800`. This fixture contains two SNVs, a tandem insertion/duplication hypothesis, and a deletion in a disconnected block. Expected: 44 fragments, four alleles, two blocks, three groups in the first block, and four collapsed group rows. The `noisy.*` fixtures contain 200 distinct SNP patterns; the strict display separates overlapping contradictions while collapsing only mutually compatible fragments. Every record is retained. The additional `paired.*` fixtures contain 12 paired fragments (24 records); `synthetic:1001-1800` collapses them into one overview row. Selecting `synthetic:1051-1110` also exercises indexed retrieval of mates beyond the selected interval.

## Interpretation and current bounds

Two-bit version 0 in either byte order is supported; version 1 is rejected explicitly.

The viewer reports supported **observed local allele patterns**, not tumor clones, cellular fractions, or a guaranteed reconstructed full haplotype. Partial observations remain partial. Ambiguous fragments are not forcibly assigned. No joint alternate reads does not prove trans. Groups do not have identities across disconnected blocks. The graph records repeated fragment co-observation, and the grouping engine deliberately retains observed patterns rather than inventing a whole-gene sequence through transitive compatibility.

Support counts are fragments (query name within sample/read group for paired alignments). Duplicate filtering is configurable for amplicon assays. UMI family consensus is not implemented. Mapped mates on the same contig are retrieved by name/read group from up to 512 indexed windows, including mates outside the selected interval. No full-file scan is performed. A mate on another contig, unmapped mate, or unavailable mate remains explicitly incomplete; it is not drawn at a fabricated position. Reference comparison and phasing remain limited to the selected interval. Supplementary records contribute within the same named template; ambiguous/chimeric mappings need manual review. Quality thresholds and minimum support are exploration settings, not validated limits of detection. Missing/255 qualities are not accepted as high confidence. Reference sequence differences are evaluated directly against the selected local reference; sequence digests are not automatically verified.

SNVs and CIGAR indels are discovered automatically. Supported exact soft-clip anchors can propose local insertion/deletion hypotheses; clips without unique local anchors remain unresolved. Local event paths use affine-gap alignment, preserve original alignments, and show reference gaps around inserted sequence. Analysis and event-table coordinates are left-normalized for equivalent indels. Pileup glyphs keep the original BAM CIGAR anchors in both collapsed and expanded rows; normalization never relocates a drawn SNP or indel. Repeat-equivalent raw indel placements are counted separately for summary display, so a placement without minimum supporting fragments is not drawn as an invented consensus anchor. The tandem duplication label is a local sequence hypothesis, not a measured copy number. This release does not perform unrestricted de novo assembly, distant rearrangement reconstruction, clonal phylogeny, or probabilistic full-mixture inference. Multi-allelic/overlapping variants remain separately flagged for inspection.

Limits: 200 kb selected interval, 100,000 returned alignment records, two million sequenced bases, 500 candidate alleles, 32 supported patterns per block, 60-second job watchdog, 64 MiB per file-range/index read, 128 MiB decompressed per BGZF query, and 256 MiB cumulative BAM I/O per query. Exceeding a limit fails the analysis explicitly. These bounds do not constitute a hard browser heap limit: parser allocations and concurrent reads can peak higher. Dense windows should be narrowed. Event reconstruction supports up to 1 kb insertion or 5 kb deletion, subject to a two-million-cell alignment budget. Event realignment processes every overlapping record with usable sequence, or fails explicitly if its cumulative 20-million-cell budget is exceeded; it never silently takes a subset. The default **Haplotype groups** view displays one collapsed consensus row per supported local group, plus compact ambiguous, uninformative, unresolved, or excluded evidence rows when needed. Compatible partial reads and endpoint or clipping differences can share a row, but any explicit A/C/G/T disagreement at an overlapping genomic position is an absolute veto, including reference-versus-alternate calls, singleton variants below discovery support, and low-quality visible disagreements. Those observations remain separate rather than being voted away. A partial observation must match a unique supported group with at least one shared informative site to be assigned; otherwise it stays in a collapsed uncertainty bucket. With no supported variants, connected covered regions are partitioned into mutually SNP-compatible coverage rows. Supported distinct structural links can subdivide a local group; weak structural annotations remain available in the structural panel without creating thousands of summary rows.

Group rows are local to a phase block; overlapping complex-allele blocks keep overlapping display scopes so their defining sites remain visible. A fragment spanning disconnected blocks can appear in each block's summary; global retained-fragment counts are unique, and per-row counts must not be summed as clone counts. Summary glyphs show unanimously observed callable alleles; overlapping mates vote once. A conflicting observation is never overridden by a majority. A discordant mate pair is isolated, and a quality-passing conflict is marked with an amber cross instead of a consensus allele. Variant discovery support thresholds remain independent from this stricter display rule. Reference-like summary shading is a coverage overview, not a claim that unobserved positions were sequenced or that an ambiguous mixture has one haplotype.

Click a summary row or its accessible table entry to expand all member pairs. Global expand/collapse controls are also available. The **Exact-pattern families** inspection mode retains the previous stricter grouping. Duplicate-marked fragments are retained in display counts and can join a compatible group, while the analysis duplicate filter remains independent. Use the scrollbar, mouse wheel, or Up/Down/Page Up/Page Down/Home/End keys to reach every row. Only the viewport is painted; no reads are downsampled. Gray reference-matching bases, colored mismatches, purple insertions, deletion strokes, and neutral clips remain distinct.

## Variants of interest

Open **Variants of interest** above the pileup and enter one variant per line, then select **Check and highlight variants**. The same list is applied when analyzing a region. The check reruns local analysis, retains requested alleles even below automatic discovery support, highlights their original positions with V1/V2 guides, and shows a compact clickable pair matrix with exact fragment counts. Click a resolved entry to center the view. Invalid lines have separate errors and do not discard other valid entries. Up to 32 entries are supported, and all must lie in the analyzed interval to be evaluated.

Supported examples:

```text
chr17:7579472G>C
chr17:g.7579472G>C
chr17 7579472 G C
NM_000546.6:c.215C>G
c.215C>G
```

Genomic substitutions use one-based positions. HGVS DNA substitutions, deletions, literal-sequence insertions, tandem duplications, and delins are supported, including coding UTR and intronic positions. Bare `c.` entries require a selected transcript; gene-qualified entries require an unambiguous transcript or an explicit selection. Versioned transcript accessions are checked against the selected assembly's embedded RefSeq catalog. CDS bounds, exon order, negative-strand alleles, and intronic anchor boundaries are accounted for. Reference alleles are validated against the local genome. RefSeq-to-genome sequence differences, unknown reference accessions, uncertain positions, protein-only descriptions, and non-literal complex insertions are not silently inferred. Ranges crossing splice junctions also require an explicit genomic DNA event. For these, use an explicit genomic allele. This is a bounded local notation resolver, not a general HGVS validation service.

**Cis supported** means alternate alleles are directly observed together on at least the configured minimum supporting fragments. **Possible cis** has direct support below that threshold. **Mixed patterns** contains joint alternates and additional separate patterns. **Unphased** means the current evidence cannot establish the relationship. Separate alternate patterns do not establish trans in a mixed tumor; an optional explicit diploid assumption permits a qualified **Trans supported** label when both reciprocal single-alternate patterns meet support and no joint-alternate observation exists. That label is conditional on the assumption and finite sampled evidence. The matrix provides both-ALT, left-only, right-only, both-REF, other-allele, unknown and jointly informative fragment counts. Paired ends count once; there is no imputation across coverage gaps.

## Exploratory group hierarchy and fractions

A compact dendrogram sits alongside the pileup on wide screens and above it on smaller screens. Block chips show each region's supported group count; select a block to compare its groups. The similarity cut slider displays an exploratory count of similarity clusters without changing phasing, display-family membership, or read counts. Hover/focus a tree group or hover a pileup group for counts, pairwise SNP comparisons, and observed site allele fractions. Click a tree leaf to locate its pileup group. Detailed phase/variant panels are collapsed by default.

The hierarchy uses [average linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) on the proportion of differing jointly observed SNP sites. Multiallelic coordinates count once. Missing calls do not count as reference matches; merges require a defined comparison for every pair across the proposed clusters. Groups without shared observed SNPs remain in separate trees. Only SNPs contribute to branch distances: indels and structural differences remain represented in groups, and groups identical at SNPs can join at zero SNP distance. The tree expresses similarity, not ancestry, mutation order, or a validated cellular clone count. Independent blocks retain their own trees and counts, so a region with seven groups and one with three groups are both represented without forcing identities between them.

**Group share** is uniquely assigned eligible fragments divided by all eligible fragments informative at any site in that block, including ambiguous/conflicting observations in the denominator. Assigned alignment-record counts are also shown. **Site VAF** is the count of eligible fragments carrying the displayed allele divided by fragments with a callable observation at that site. Paired ends count once. A variant shared by two groups can have a higher VAF than either group's fragment share. Ambiguous fragments are not allocated to a group, excluded/duplicate-filtered records do not contribute to these fractions, and neither statistic estimates tumor-cell fraction. Hover details expose every numerator and denominator.

## Split alignments and rearrangement evidence

Using the [SAM tag definitions](https://samtools.github.io/hts-specs/SAMtags.pdf), the structural panel inspects local records' `SA` tags, available supplementary alignments for the same read segment, and mapped-mate coordinates/orientations. It reports possible interchromosomal translocations, same-chromosome split links, inversion-like complementary split alignments, same-strand mate clues, and distant same-chromosome mates. Signals require an anchor overlapping the selected locus. Remote regions are not fetched or assembled. A gene lookup button uses the selected embedded assembly catalog to annotate genes overlapped by the partner alignment; overlap alone does not establish a fusion.

Split geometry accounts for strand, soft/hard clipping, and query-segment overlap. Overlapping alternative mappings are labeled ambiguous rather than called rearrangements. Nearby anchors are clustered within 100 bp of a seed, with full coordinate ranges retained. More-than-10-kb mate separation is a display threshold, not an empirically fitted insert-size cutoff. Same-strand mate clues are library dependent; no library orientation is assumed or calibrated. Mate-only coordinates are alignment starts, not resolved breakpoints. Missing `SA` tags or mate metadata limits what can be identified; absence of a signal cannot exclude a structural event.

Counts distinguish total supporting fragments, fragments passing both local/partner mapping-quality filters, unknown partner MAPQ, and duplicate-marked fragments. Mate MAPQ is taken from `MQ` or an already available mate record. Primary/supplementary records and the two ends of a fragment do not get independent support votes. Singletons and low-quality signals remain visible with their counts. Families with different structural links are kept separate and marked **SV** in the pileup.

The `structural.*` synthetic fixtures contain split-read translocation, inversion-like split alignment, and cross-chromosome mate examples with three fragments per signal. Select `synthetic:1001-1800`. These annotations are local supporting evidence, not genome-wide structural variant calls or reconstructed partner sequences.

## Privacy

All selected-file access and analysis happen in the browser. There are no uploads, runtime downloads, remote references, telemetry, saved preferences, automatic exports, or session persistence APIs. Scripts and styles are hash-authorized by CSP; connections, frames, forms, objects, and external assets are blocked. The dedicated Blob worker is terminated on reset/cancellation. Worker/network and actual standalone launch paths are included in browser checks. Read names are exposed only through explicit record inspection. Public reference catalogs are embedded; reference bases are selected locally every session.

This describes application behavior, not forensic erasure. Browsers/OSs control memory reclamation, swap, extensions, crashes, and recovery. Use a trusted browser/profile and appropriate device controls for sensitive data. The application is not clinically validated.

## Verification

Verified on November 10, 2026 in desktop Chrome 157.0.8103.50 on macOS. The code and browser regression suites cover analysis, paired display families, collapse/expand, and full vertical scrolling. The standalone browser regression passes synthetic demo analysis, multi-block BAM decoding with BAI/FASTA and CSI/two-bit, gene lookup, expanded event alignment, clear/reload, and checks for page errors and outbound requests. No real patient BAMs were used. Firefox, clinical accuracy, high-depth performance targets, and browser/OS forensic erasure have not been validated.

## Development

```sh
npm ci
npm test
npm run fixtures
npm run build
npm run test:browser
```

Only development requires Node/npm. The browser test uses the installed Chrome binary on macOS; set `CHROME_PATH` for another installation. It exercises the actual `file://` artifact using synthetic BAM/BAI and FASTA/FAI, verifies event alignment and gene lookup, and checks page errors, outbound requests, reset, and reload.

`src/analysis.js` contains reference projection, indel normalization, fragment evidence, conservative pattern grouping, anchored clip hypotheses, and affine alignment. `src/worker.js` owns BAM/reference access and jobs. `src/io.js` exposes local File reads only. `src/app.js` owns controls and Canvas rendering. `scripts/build.mjs` replaces the upstream optional WASM decompressor with `src/bgzf.js`, bundles everything, and generates CSP hashes and `dist/SHA256SUMS`. The release contains no WebAssembly. Build inputs are recorded in `dist/build-inputs.json`.

Public hg19/hg38 gene/transcript catalogs are compressed in `references/*.js`. Their source URLs, retrieval date, and source checksums are in `references/manifest.json`. They are NCBI RefSeq coordinates distributed through UCSC. No case data belongs in this repository. Dependency attribution is in `dist/THIRD_PARTY_NOTICES.txt`.

See `IMPLEMENTATION_PLAN.md` for the broader design, including extensions beyond this release's conservative observed-pattern reconstruction. Actual patient BAM testing is intentionally left to the user.
