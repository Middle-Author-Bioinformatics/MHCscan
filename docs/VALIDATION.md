# Validation and scientific limits

Validation was performed in a Linux environment using Python 3.12, BioPython 1.88,
PyHMMER 0.12.3, Pyrodigal 3.7.1 and pytest 9.2.0. This is functional and logic
validation, not a measured sensitivity/specificity benchmark for Gram-positive EET.

## Automated checks

Run core tests:

```bash
python -m pip install '.[analysis,test]'
python -m pytest -q
```

Run the externally sourced integrations with the documented checkouts and database:

```bash
MHC_EET_FEGENIE_REPO="$PWD/external/FeGenie" \
MHC_EET_FEET_REPO="$PWD/external/FEET" \
MHC_EET_SURFACE_DB="$PWD/databases/surface.hmm" \
python -m pytest -q
```

The full local suite passed **24 tests** (21 core/cache tests and 3 external
integrations). Without the three optional environment variables, those external
integration tests are explicitly skipped. Core tests do not fetch scientific data.
A GitHub Actions workflow is included for Python 3.10/3.12; the remote workflow's
result is not claimed as part of local validation.

Covered behavior:

- One/two/five canonical motif thresholds; mixed canonical/noncanonical motifs;
  MHCscan legacy counts; FEET's literal long-regex behavior and the alternative
  negation interpretation; isolated/internal LPXTG versus terminal architecture.
- Exact protein/CDS mapping, Parent mappings, ambiguity errors, shuffled GFF/order,
  GenBank translations and partial CDSs, all-CDS counting including untranslated
  pseudogenes, no cross-contig associations and no suffix-based ID assumptions.
- Eight CDSs on each side, immediate/intervening counts, signed gap and overlap,
  strand-relative context, convergent/divergent orientation, gene/bp-window
  intersection and bp-only selection.
- Linear edge truncation, explicit circular wrap and boundary spacing, partial
  anchors, unavailable coordinates and conservative disabling around compound CDSs.
- HMMER repeated hits, domain coordinates/coverage, GA thresholds, missing-GA errors,
  explicit E-value fallback, cached-result reuse and corrupted-raw invalidation.
- Same-protein domains versus separate sortase machinery; MHCs without qualifying
  neighbors retained; missing localization and extracellular exposure remain
  unavailable; source-preserving FeGenie/FEET/MHCscan imports.
- Native Phobius, PSORTb, SignalP and CELLO Results parsing; unknown imported protein
  IDs rejected; independent localization feature records retained.
- Deterministic Pyrodigal prediction and mapping from a synthetic nucleotide input;
  prevention of silent mixing with existing annotations.
- Sequence-only CLI, no-candidate output, completed-run resume, corrupted output
  detection, offline HTML tables/embedded downloads/SVG presence.

## External analyses actually run

1. Downloaded all six default HMMs directly from official InterPro endpoints and
   verified their current metadata/headers, accession versions and GA cutoffs.
   Official Pfam release metadata reported 39.0 / 2026-09 / UniProtKB 2026_02.
   Tested extraction of the saved HMM collection without network access.
2. Searched model-derived two-repeat positive-control sequences against all six
   official Pfam models with GA thresholds. Each model returned at least two own-
   family hits. These controls test repeated-hit handling; they are **not real
   proteins or a biological validation set**.
3. Ran all **43 FeGenie category/model entries** in iron_reduction and iron_oxidation
   against the three Thermincola reference proteins in FeGenie's
   `TherJR_SLCs.faa`, and ran the six Pfam searches on those same proteins.
   FeGenie returned family evidence for two references. No default surface-domain
   hits were returned at the supplied GA thresholds; no domain hit was invented
   to match an expected localization. There were 49 model-search entries total.
4. Executed the original FEET `Find_multihemes_CNO.py` with parameters `3 5`, from
   the inspected external checkout, on both synthetic motifs and real references.
   Literal code-compatible counts and pass/fail agreed, including long-motif
   discrepancy fixtures. Original count/list/FASTA outputs were retained.
5. Generated a complete synthetic matched-annotation report with repeated
   same-protein SLH hits, separate sortase and G5 neighbors, unrelated CDS context,
   and a partial MHC at a fragment edge. The imported hits are explicitly marked
   `synthetic_fixture_NOT_real_search`, never presented as real domain predictions.

### Thermincola reference interpretation

| Reference FASTA ID | Original annotation | Canonical count verified |
|---|---|---|
| 646799199 | YP_003641333 / TherJR_2595 | 9 |
| 646797728 | YP_003639887 / TherJR_1122 | 6 |
| 646796949 | YP_003639120 / TherJR_0333 | 10 |

TherJR_2595 corresponds to OcwA, a characterized surface MHC with EET-related
biochemistry. The other proteins are pathway-associated/proposed wall/periplasm
components; they do not all have equivalent experimental evidence. See
[*Crossing the Wall*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7911101/) for the
Thermincola pathway and characterization context. These proteins were tested as
sequence controls, not as a claim that motif counts establish their locations.
No matching coordinates were supplied in the sequence-control run, so their
neighborhoods were correctly unavailable. No localization predictions were supplied;
localization was correctly unavailable despite literature context.

Reproduce that sequence-control run:

```bash
mkdir -p /tmp/thermincola_eet_reference
cp external/FeGenie/hmms/iron/iron_reduction/non-aligned/TherJR_SLCs.faa \
  /tmp/thermincola_eet_reference/Thermincola_reference.faa
mhc-eet run --input /tmp/thermincola_eet_reference --input-type proteins \
  --surface-db databases/surface.hmm --fegenie-repo external/FeGenie \
  --feet-repo external/FEET --output thermincola_eet_reference --threads 2
```

Generate the synthetic annotation/report fixture:

```bash
python examples/make_demo.py /tmp/mhc_demo
mhc-eet run --input /tmp/mhc_demo/proteins --input-type proteins \
  --evidence-dir /tmp/mhc_demo/imports --skip-surface \
  --output /tmp/mhc_demo/results
```

## Not run / limitations

- Full FeGenie clustering/DIAMOND validation and the complete legacy MHCscan
  Phobius/HHsuite workflow were not run; the wrapper specifically uses their
  documented subset and supports imports of existing outputs.
- Live CELLO submissions, PSORTb inference, SignalP inference, lipoprotein or
  separate TM topology prediction were not run. Parsers were tested with explicit
  fixtures; missing predictions are unavailable. Licensing/access and heavyweight
  predictor models are intentionally outside the core dependencies.
- The original FEET manual CELLO/porin pipeline was not run or required, and no
  Gram-negative eligibility rule was used to exclude Gram-positive candidates.
- The offline report's generation, embedded downloads, tables and SVG structures
  were checked. A headless Chromium interactive check was attempted but blocked
  by the execution environment's socket restriction; browser-rendered appearance
  and interactive search were not claimed as visually validated.
- No complete experimental Gram-positive positive/negative benchmark was available
  in the inspected repositories. No real whole-genome cohort performance or
  macOS execution was benchmarked. Genome calling was checked using synthetic
  nucleotide input; existing biological references tested sequence analysis.
- Compound-CDS contigs retain sequence analysis but disable neighborhood inference.
  Noncanonical sortase signals, distant accessory loci and non-MHC EET can be missed.

All resulting scientific annotations are candidates or source-specific predictions.
None is an experimental demonstration of EET, extracellular exposure, interaction,
or an operon.
