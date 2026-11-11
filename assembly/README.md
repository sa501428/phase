# Phase assembly

Open [index.html](index.html) directly in desktop Chrome or Edge for the assembly-first viewer. The automatic grouping variant is available at [../automatic/index.html](../automatic/index.html).

The running viewer was updated on November 11, 2026 with configurable **Minimum variant VAF (%)**, default **2**, in the assembly/evidence settings. Alleles below the threshold are hidden and excluded from matching, graph variant evidence and clash checks. Exactly 2% is included; 0 disables suppression. Changing the cutoff automatically refreshes the loaded reads and variant displays. The original reads remain available. Without a reference, the modal observed base supplies the baseline.

The original `source.tar.gz`, documentation, literature notes and notices remain unchanged. That historical source predates the VAF update. Maintained source is in `v2/` and shared `src/` modules in the current repository; run `node scripts/build.mjs v2` to rebuild this HTML and update its checksum. The checksum manifest verifies the maintained HTML alongside the original snapshot files. No private case files are included.
