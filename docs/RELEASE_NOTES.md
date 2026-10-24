# Release preparation: 2026-10-24

The archived numerical results are from the completed 30-case species-dependent
external-transfer screen, plus six verification cases. No scientific equation,
physical coefficient, FE value or saved concentration/isotope array was changed
for this release.

Changes from the working calculation files:

1. Restored the required `model/` and sibling `fe_audit/` input layout, making
   solver execution self-contained.
2. Kept only numerical FE fields; internal workbook paths/cell annotations are
   removed. The numerical input values are unchanged.
3. Changed new-run hash keys from absolute workstation paths to repository-relative
   names; sanitized archived JSON path keys similarly, retaining the input hashes.
4. Added portable screen/plot/check/post-processing commands. The older
   Windows-only runner, logs, locks and HTML-integration side effects are omitted.
5. Preserved the existing Figure 4 redraw, with PDF output opt-in rather than
   automatic. No figure-layout redesign or experimental-image editing was done.
6. Updated the parameter ledger's old base-module numerical defaults to the
   actual corrected solver (48/96/24 grid, BDF 1e-9/2e-10 and tightened closure).
7. Added explicit warnings about effective versus actual flow, pH input 3.4,
   approximate concentrated-electrolyte chemistry, restricted convergence checks,
   and depleted-liquid isotope precision.

`provenance/imported_files.json` records original and public file hashes.
`provenance/original_input_hashes.json` identifies the four original campaign
inputs. A public checksum manifest covers the release files (excluding generated
run outputs and itself). Portability additions do not retroactively claim the
whole screen was newly recalculated on the release date.

This release contains no manuscript/SI documents, experimental raw archive,
literature PDFs, author emails, comments, correspondence, local HTML reports,
credentials or Git history from the larger workspace. NMR/GC–MS processing and
legacy active-learning notebook repairs remain outside this spatial release.
