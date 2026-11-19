# Executed validation

- Ten Node core tests pass, including exact source-file hashes, derived coordinate agreement, proper rigid transformations versus reflections, malformed/sparse inputs, 214 matched residues, declared domain partitions, and alignment-invariant internal distances.
- Independent NumPy 2.2.6 Kabsch SVD checks both real fits and 30 seeded rigid-transform cases. Maximum coordinate discrepancy is below 1e-8 angstrom in those fixtures; exact results are in numpy-comparison.json. SciPy 1.15.3 is present in the optional recorded environment but is not the alignment oracle; NumPy SVD is.
- The full browser workflow passed locally and against the public deployment, using an isolated Chromium headless shell with software WebGL. Desktop 1440x1000 and an emulated 390x844 phone are covered. Tests execute endpoint labels, interpolation disclosure, alignment sensitivity, region selection, atom/backbone mode, overlay, ruler, JSON and PNG downloads, reduced motion and a forced WebGL-failure path. See browser-results.json.
- The first phone check failed: the hero canvas exceeded a 390px viewport by 7.5px. Diagnosis identified the oversized cover-specimen rule. Reducing its width and offset removed the overflow; both later local and production runs passed the same check. This was a layout failure, not a detector or structural-alignment failure.
- Archive integrity and extraction are checked; the extracted zero-dependency CLI reproduces the endpoint ruler. Dependencies, environment files and deployment credentials are excluded from source downloads.

Not executed: physical iPhone/Android tests, native app tests, Safari/Firefox coverage, clinical/experimental validation, molecular dynamics or mechanistic inference. Numerical agreement is not evidence that interpolated frames are physically realized. The original structural determinations belong to their depositors, not this project.

## Version 1.1 protein workspace

All 26 core tests, the original browser regression suite, and the expanded protein-workspace browser workflow passed locally. Live RCSB CORS downloads/search were exercised, not simulated. See protein-workspace.md and protein-browser-results.json. The new construction tests validate declared geometry and template provenance; they do not validate a native fold or physical plausibility.
