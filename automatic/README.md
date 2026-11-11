# Phase automatic grouping

Open [index.html](index.html) directly in a browser for the pre-assembly viewer. It is self-contained, with no server or downloads required.

The running viewer was updated on November 11, 2026 with a configurable **Minimum variant VAF (%)** under **Evidence filters**, default **2**. Alleles below the threshold are hidden and excluded from matching and clash checks. Exactly 2% is included; 0 disables suppression. Changing the cutoff automatically refreshes the loaded reads and variant displays. The original reads remain available.

The original source snapshot, documentation and notices remain unchanged in [v1-2026-11-10](v1-2026-11-10/). That historical source predates the VAF update. The maintained source is `src/` in the current repository; run `node scripts/build.mjs src` to rebuild this HTML and update its checksum. The checksum manifest verifies the maintained HTML alongside the original snapshot files. No patient files are included.
