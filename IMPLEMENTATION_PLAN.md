# Phase manual haplotype workspace

Implemented at top-level `index.html`. The assembly-first viewer and its frozen source are preserved at `assembly/index.html` / `assembly/source.tar.gz`; the original viewer remains at `automatic/index.html`.

- Explicit gene/transcript/exon/coordinate target; indexed BAM/BAI/CSI plus local 2bit or FASTA/FAI reference.
- IGV-style reference/coverage tracks and densely packed paired alignments with independent item hit testing; colored SNPs/indels, visible soft clips, MAPQ/CIGAR/strand details through right-click with optional hover off by default, zoom, horizontal drag/pan and vertical row scrolling.
- Whole-pair click selection and Shift-click ranges across scrolling; named manual merges, expandable original members, ungroup/dissolve, reversible hide/restore, rename and Undo/Redo.
- Contradictory callable bases/deletions/insertions veto merges. Unknown gaps remain unknown; quality changes flag conflicts without rewriting assignments. No automatic grouping or read downsampling.
- Average-linkage SNP similarity for constructed haplotypes, transparent group pair shares, site VAFs and within-group allele fractions. Hidden-evidence behavior is explicit and configurable.
- Validated assignment JSON save/restore keyed to the loaded evidence; observed-block FASTA export with no reference-filled gaps and no arbitrary conflict consensus.
- Synthetic model and browser regressions, including ten unequal nested/sibling somatic lineages, paired MAPQ/soft clips, strict conflicts, missing SEQ, hide/restore, history and local scope.
- Self-contained file:// release, embedded offline annotation, RAM-only state, CSP network blocking and explicit resource limits.

The workflow supports human construction of local sequence hypotheses. It does not claim a statistically inferred number of cancer-cell clones, cellular abundance, tumor ancestry or remote structural-event reconstruction.

Packed-layout updates: nonoverlapping paired templates reuse lanes, larger configurable glyphs, duplicate-flag filtering by default, Shift-drag area selection, independent collapsible side panels, compact top controls and a viewport-filling pileup. VAF, consensus and coverage respect the explicit duplicate toggle.
