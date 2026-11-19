# Protein workspace — 1.1

The original adenylate-kinase comparison remains available. `explore.html` is a separate general coordinate viewer and idealized sequence constructor; it does not apply the original study's CORE/NMP/LID domain boundaries to other proteins.

## Sources and parsing

RCSB protein-name search uses the official Search API's full-text service, limited to six experimental-entry results; an entry's metadata is read through the Data API. A classic or extended PDB identifier downloads uncompressed mmCIF through the official file service. Browser requests are restricted to files.rcsb.org, search.rcsb.org and data.rcsb.org. Searches/downloads are explicit user actions; sequence construction and local coordinate parsing do not send sequence/file contents to these services.

The importer uses first model/data block, author chain/residue labels with label fallbacks, insertion codes, and blank/A alternate locations with blank preferred. It retains protein heavy atoms, known modified amino acids and otherwise recognizable ATOM peptide backbones. Water, ligands and nucleic acids are omitted. This is not a full crystallographic reader, assembly generator, chemical validator, confidence evaluator or sequence aligner. BinaryCIF, compressed files and arbitrary archive variants are not supported.

Limits are 12 MiB input, 100,000 parsed protein atoms and 24,000 displayed atoms. If an entry exceeds the display limit, a supported single chain is selected and identified. No model is silently truncated to the atom limit. The original coordinate bytes remain downloadable, including data excluded from the view. Hashes are computed on those original bytes. Coordinate-derived FASTA omits missing residues and labels that limitation.

PDB and mmCIF encodings of 1CRN agree atom by atom in the retained fixture: 327 heavy atoms, 46 Cα residues. The 4HHB mmCIF fixture retains 4,384 protein atoms and 574 Cα residues across four chains; heme and water are absent from the display but remain in an original-file download. Fixture sources: https://files.rcsb.org/download/1CRN.cif, https://files.rcsb.org/download/1CRN.pdb, https://files.rcsb.org/download/4HHB.cif; retrieved 2026-11-18.

## Sequence construction

One FASTA record or 3–500 letters from the 20 standard amino-acid alphabet is accepted. Coordinates are constructed from declared bond lengths/angles and chosen uniform φ/ψ torsions, with ω fixed at 180°. Presets use α-helical (-60°, -45°) or extended-strand (-135°, 135°) torsions. The backbone is centered before export so the advertised maximum extended chain fits legacy PDB coordinate fields.

Optional side chains use 20 bundled Chemical Component Dictionary ideal-coordinate templates. They are rigidly attached through their N/CA/C anchors, using the existing verified Horn fitter. A single reference side-chain shape is used, not rotamer prediction or energy optimization. Template source URLs, retrieval times, hashes and anchor-fit RMSDs are retained. Raw component CIF fixtures are included for independent regeneration; the current official CIF definitions include the ideal-coordinate columns. Legacy `_ideal.pdb` retrieval was unavailable during preparation and is not used.

No sequence-specific native fold is inferred. No hydrogens or terminal caps are added. Steric clashes, ring closure, protonation, side-chain/backbone junction validity and residue-specific torsion preferences are not validated or optimized. Proline/glycine-specific restrictions are not enforced. PDB exports carry explicit construction/non-prediction remarks. This is a geometry-building feature, not a substitute for folding, force fields, molecular dynamics or structural validation.

## Executed checks

All 26 Node tests pass: the original comparison tests plus real multichain imports, independent PDB/mmCIF coordinate agreement, malformed input rejection, bounded stream/display behavior, FASTA validation, backbone bond lengths and torsions, source hashes for all 20 CCD templates, side-chain identity/internal-distance preservation, and 500-residue PDB round trips.

The protein browser workflow exercises real RCSB download/search, local file imports, original-file preservation after chain filtering, sequence construction without external requests, PDB/JSON/FASTA-related state, invalid-input recovery, extended torsions, desktop/390px layout and a forced non-WebGL projection fallback. Exact executed checks and environment are in protein-browser-results.json. The original comparison browser suite also passes.

Browser checks use software WebGL and an emulated phone viewport, not physical-device validation. Bone/atom geometry checks are numerical, not proof of molecular plausibility. The sequence constructor is consistently labeled as constructed geometry, not a predicted fold.

Primary documentation: https://www.rcsb.org/docs/programmatic-access/file-download-services, https://search.rcsb.org/, https://data.rcsb.org/, https://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html.
