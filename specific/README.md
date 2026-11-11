# Phase — specific variants

Open **[index.html](index.html)** directly in desktop Chrome or Edge. This release is a single self-contained HTML file containing all JavaScript, CSS, BAM reading libraries and offline transcript annotations. No server, installation, reference-genome file or network connection is needed to use it.

1. Choose a coordinate-sorted BAM and its matching BAI or CSI index, then **Open files**.
2. Choose a single sample if the BAM contains several. Samples and read groups remain separate; read pairs are never assembled across read groups.
3. Paste variants into the large text box and **Check variant list**. Correct flagged lines. Repeated minimally normalized alleles are coalesced. There is no fixed variant-count limit; larger lists and results are paginated.
4. Set quality and optional helper-SNP thresholds, then **Phase specific variants**.
5. Inspect a requested pair to see its frequency denominators, per-site observations, helper path if available, and the contributing R1/R2 pairs. Switch the table to **All requested pairs** to include pairs without direct joint coverage.

## Paired-end frequency counts

Each primary physical read pair contributes once, including when one end covers A and the other covers B across a long unsequenced insert. An indexed query at each requested site retrieves the relevant ends; helper mode also queries bounded neighborhoods and intervals between nearby targets. There is no BAM-wide scan. Overlapping mates count once. High-quality mate disagreements become conflicts and are excluded from the frequency denominator. Inconsistent primary-end identities or reciprocal mate coordinates fail explicitly. Cross-contig discordant pairs do not establish local cis phase.

For each requested pair, the app reports:

- **Both ALT**: A ALT and B ALT on the same fragment.
- **A only / B only**: one specified ALT and the other specified REF.
- **Neither**: both specified REF.
- **Other / conflict**: callable third alleles or disagreement within the fragment, excluded from REF/ALT frequencies.
- **Together frequency**: both / (both + A-only + B-only).
- **Separate frequency**: (A-only + B-only) / (both + A-only + B-only).
- **Conditional fractions**: B ALT among A ALT, and A ALT among B ALT.

Every percentage shows its numerator and denominator. Neither-ALT observations remain visible but do not dilute the ALT-bearing together/separate denominator. No coverage, low-quality bases, splice skips, incomplete indel flanks and absent sequence stay unknown. A zero denominator is shown as unknown, not zero percent. Per-site callable ALT fractions are also shown in the inspector. Separate patterns are not a claim of diploid trans phase, cellular abundance or tumor ancestry.

The paired view keeps both ends in a single evidence row, with R1/R2 coordinates, CIGARs, mapping qualities and strands. Bars share a genomic axis and the insert gap remains a dashed connector. Requested sites use teal; helper sites use purple.

The separate **Minimum variant VAF** setting defaults to 2% for supplied variants and other phase markers. It automatically removes alleles below the cutoff from rendered results, pair matching and helper links, retains exactly 2%, and can be set to 0 to disable suppression. This per-site fraction uses ALT / (ALT + REF + other), excluding conflicts and unknown coverage. The completion status reports how many variants were suppressed.

## Optional helper SNPs

Helpers extend marker relationships while retaining the requested variants as the primary output. Defaults are:

| Setting | Default |
| --- | --- |
| Minimum unique molecules carrying the helper allele | 5 |
| Minimum helper allele fraction | Strictly greater than 2% |
| Maximum gap between target neighborhoods to search | 50,000 bp |
| Minimum per-link phase-pattern agreement | 95% |

All four settings are adjustable; molecule support has a hard minimum of five. Both observed alleles need the molecule-support threshold. Helpers are biallelic SNP markers; positions in or immediately adjacent to requested events are excluded. The search includes 250 bp around targets and fills gaps up to the selected maximum, rather than searching the whole BAM.

**Molecule counting:** Read names are scoped by sample and read group. Helper discovery and links additionally collapse identical paired alignment coordinates, with RX (or UB) tags distinguishing molecules when present. Single-end records use alignment start/end/strand. Disagreements within a collapsed molecule are excluded. A duplicate flag on either end removes the whole template from helper evidence, even when duplicates are enabled for direct counts. This is a conservative coordinate/UMI molecule proxy, not UMI consensus or proof of physical independence; identical-endpoint untagged molecules can be collapsed. Distinct UMIs are accepted as supplied, without error correction.

**Reference-free allele fraction:** Without a reference genome, the modal observed A/C/G/T base is the baseline and the nonmodal base is the helper allele. The reported helper fraction is nonmodal-molecule count / callable-molecule count. It can differ from genomic-reference VAF, particularly if the genomic ALT is the majority allele. Markers with more than two observed base alleles are not used.

**Links:** A link must observe both phase patterns with at least the minimum molecule support: ALT/ALT plus REF/REF for a coupled relationship, or ALT/REF plus REF/ALT for an opposite relationship. The supported orientation must meet the agreement threshold. Links with third alleles or molecule conflicts are excluded. Edges must include a helper SNP; the app does not silently chain only the supplied variants.

Purple markers in an assisted path are explicitly named, with per-link molecule counts and agreement. Click a link to inspect its contributing original read-pair templates. Those rows can outnumber the coordinate/UMI-deduplicated molecule counts. Contradictory signed cycles make the reachable helper block unresolved. The app never multiplies local percentages into an invented end-to-end frequency; a helper chain is indirect phase evidence, with its own limitations in mixed samples.

## Accepted input

Input positions are **1-based**. The internal implementation uses zero-based coordinates and insertion boundaries.

```text
chr17:7674220C>T
chr17:7674220:C:T
17-7674220-C-T
chr17 7674220 C T
chr17:g.7674300_7674302del
chr17:g.7674300_7674301insAC
NM_000546.6:c.743G>A
```

Supported tables include VCF (with or without its column header), four-column CHR/POS/REF/ALT tables, and headed CSV/TSV/MAF with recognized chromosome, position, reference and alternate columns. MAF uses `Chromosome`, `Start_Position`, `Reference_Allele`, and `Tumor_Seq_Allele2`. Multi-allelic ALT cells are split; quoted CSV cells are supported. A headed table should be pasted as a single table. `-` represents an explicitly empty allele at a 1-based event coordinate; VCF-anchored indels are recommended.

HGVS supports substitutions, del, ins, dup and delins. Transcript mapping uses the embedded GRCh37/hg19 or GRCh38/hg38 RefSeq catalog; choose the build matching the BAM. Exact accession versions are respected, ambiguous gene/transcript mappings fail, negative-strand alleles are reverse-complemented, and intronic positions are mapped only at valid exon boundaries. Coding events spanning separate exons are rejected as ambiguous genomic events. Genomic `NC_...version` accessions are checked against the selected primary-chromosome assembly mapping. Protein-only annotations, rs IDs, structural/symbolic alleles and uncertain coordinates require explicit genomic alleles.

No independent genome is available to validate the supplied REF or assembly. **Indels are assessed at their stated alignment coordinates** after trimming shared allele anchors; equivalent shifted CIGAR events in repeats may be missed. Sequence-free HGVS deletions use intact span versus exact deletion/replacement evidence. Sequence-free duplications require the copied source and inserted sequence to be observed on the same alignment. These are explicit limitations of BAM-only operation.

## Local processing and limits

All case data is held in browser memory. There are no uploads, application servers, runtime downloads, cookies, localStorage, sessionStorage, IndexedDB, service workers, saved sessions or result files. Clear session terminates the worker and removes case inputs and results; leaving/reloading the page clears the app state. A strict Content Security Policy blocks runtime network connections. The build scripts and regression checks are developer tools, not a runtime server.

Defaults exclude duplicate-flagged templates, secondary/supplementary/QC-failed/unmapped records, MAPQ 255 and missing/255 base quality. MAPQ and base-quality thresholds default to 20. Missing or malformed read sequence does not vote. The duplicate toggle affects direct counts only.

Resource bounds fail explicitly, with no sampling or partial final results: 100,000 fetched alignments, 20 million read bases, 100 kb per read, 10 kb per requested allele, 256 MiB cumulative BAM local I/O, 64 MiB per index/text file, 128 MiB per BGZF decompression operation, two million searched helper reference bases, one million helper-covered positions, 2,000 qualifying helper SNPs, two million observed site-pair cells, 20 million fragment/site-pair comparisons and a three-minute inactivity watchdog. BAI positions must be within its 2^29 coordinate capacity; CSI is supported for larger coordinates. The pair table enumerates unobserved requested pairs lazily rather than allocating a full quadratic matrix. No processing results are persisted.

## Development and validation

From the repository root:

```sh
npm run build:specific
node --test tests/specific.test.js tests/vaf-specific.test.js
npm run test:browser:specific
```

The build writes only this version's standalone HTML, checksum and third-party notices. It does not replace the top-level viewer. Source is in this directory; shared local file, BGZF and CIGAR readers are reused from `src/`. Offline transcript mapping is isolated in `transcripts.js`. The main thread renders paginated tables; a Blob worker does file reading, variant parsing, molecule counting and helper discovery. Processing can be canceled by terminating the worker.

Synthetic tests cover formats and coding coordinates, >32 variants, BAM paired-end bridging, reference-free SNP/indel calls, quality and duplicate handling, conflicting overlaps, coordinate/UMI-deduplicated helpers, thresholding, chain provenance, contradictory cycles, bounds, embedded-CSP hashes, file:// operation, BAM/BAI/CSI, reset, mobile layout and absence of network/storage calls. They do not validate clinical performance or every real-world aligner/variant representation.

BAM/CIGAR interpretation follows the [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the supported annotation subset follows the [HGVS DNA recommendations](https://hgvs-nomenclature.org/main/recommendations/DNA/delins/). The repeat-positioning limitation is described in the [HGVS 3′ rule guidance](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/main/docs/recommendations/checklist.md). Annotation provenance is in [the reference manifest](../references/manifest.json); third-party licenses are in [THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt).
