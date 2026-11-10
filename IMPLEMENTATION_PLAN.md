# Assembly-first Phase: implemented scope

The current entry point is `index.html`; the previous working viewer is frozen at `archive/index.html`, with its source and documentation snapshot under `archive/v1-2026-11-10/`.

- Indexed BAM/BAI/CSI input, explicitly selected gene/exon/interval; optional local 2bit or FASTA/FAI reference.
- Sequence-built strand-symmetric de Bruijn graph, observed adjacency provenance, compact unitigs, folded reverse complements and GFA export.
- Bounded read-observed bubble alleles, physical-fragment linkage, minimum-support block cuts, unrestricted local pattern count, strict contradictory-SNP guard, and explicit shared/weak evidence.
- Main graph with pan/zoom, selectable local patterns, sequence differences and support/fraction inspector; paired read expansion with vertical virtualization and horizontal panning.
- Manual genomic/HGVS variant linkage, local inversion/translocation/split/mate clues, offline transcript catalog and single-exon selection.
- Synthetic regression coverage including ten unequal, nested/sibling somatic lineages, unknown/shared founder evidence, indels/repeats, weak links, missing SEQ and scope boundaries; no patient BAMs required.
- Self-contained local release, RAM-only case state, CSP network blocking and explicit limits without downsampling.

Algorithm rationale and limitations are documented in `LITERATURE_REVIEW.md`. Published methods inform the design; their complete inference engines and accuracy claims are not imported. In particular, this is an exploratory local compatible-pattern heuristic, not validated cellular clone-number inference, tumor phylogeny reconstruction, or a full remote structural-variant assembler.
