# Gram-positive candidate EET wrapper

This extension adds `mhc-eet` without replacing the existing MHCscan command.
It screens candidate multiheme cytochromes broadly, searches all proteins for
surface-related domains and FeGenie families, and describes genomic context.
It never labels these candidates as experimentally confirmed EET proteins.

## Installation

Python 3.10+ is required. Linux was tested; the dependencies also publish macOS
packages, but macOS was not tested here. Use a fresh environment:

```bash
git clone --branch codex/gram-positive-eet-wrapper https://github.com/Middle-Author-Bioinformatics/MHCscan.git
cd MHCscan
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[analysis]'

git clone https://github.com/Middle-Author-Bioinformatics/FeGenie.git external/FeGenie
git clone https://github.com/Middle-Author-Bioinformatics/FEET.git external/FEET
mhc-eet setup-db --output databases
```

No system HMMER/Prodigal executable, R, DIAMOND, BLAST, HHsuite database, CELLO,
Phobius, or localization model is required for this command. PyHMMER runs HMMER
in process; Pyrodigal runs Prodigal in process. `pip install .` installs only the
sequence/import/report core; `.[analysis]` adds both analysis dependencies.
Install the analysis extra when running HMM searches or predicting genes.

For exactly the externally inspected versions, pin the companion repositories:

```bash
git -C external/FeGenie checkout 034f554fb12528b43f7515b5c066de610a6ec6e8
git -C external/FEET checkout bd256f4df08fc70b69ede210d75f46466e5d002c
```

`requirements-tested.txt` records the Python analysis versions used for validation.
The manifest records the versions actually installed at runtime.

## Run on a folder of assembled genomes

Keep one assembly per file in an input folder (`.fna`, `.fa`, or `.fasta`):

```bash
mhc-eet run \
  --input /path/to/genomes \
  --input-type genomes \
  --output eet_results \
  --surface-db databases/surface.hmm \
  --fegenie-repo external/FeGenie \
  --feet-repo external/FEET \
  --threads 16 \
  --resume
```

Reproducible Pyrodigal **meta mode**, restricted to translation table 11, is the
default. Meta mode works with fragmented assemblies and does not require a
minimum genome length for training. `--translation-table 4` is supported explicitly;
translation tables are not inferred from organism names. Partial gene calls are
retained and flagged. Input files must be uncompressed.

This command does not infer whether an organism is Gram-positive. Select genomes
with the intended cell-envelope biology yourself. Absence of these candidates
is not evidence against flavin-mediated or other MHC-independent EET.

### Existing annotations and sequence-only mode

- `--input-type genbank`: `.gb`, `.gbk`, `.gbff`, `.genbank`. Existing CDS translations
  are used; no genes are predicted or silently translated with a guessed code.
  All CDS features, including pseudogenes without translations, enter gene order.
- `--input-type proteins`: `.faa`/`.pep` (or explicitly selected `.fa`/`.fasta`),
  with exactly one same-stem `.gff` or `.gff3`. Protein IDs must match exactly.
  `protein_id`, `ID`, `locus_tag`, and `Parent` are considered; Parent's explicit
  `protein_id` can resolve transcript IDs. Ambiguous matches fail. Use
  `--gff-id-attribute protein_id` (or another attribute) to select one mapping.
- Without a matching GFF, proteins are screened and neighborhoods are explicitly
  unavailable. Unmapped proteins in a partly mapped input are retained and flagged.
- `auto` recognizes `.faa`/`.pep` as protein input and GenBank suffixes as annotations;
  `.fa`/`.fasta` mean nucleotide assemblies. Mixed annotation types and duplicate
  genome stems fail. Adjacent annotations in genome-prediction mode cause an error;
  place the intended inputs in a separate folder or choose their annotation type.

Mapping is exported in `genomes/<genome>/cds_mapping.tsv`, including CDSs without
proteins. IDs are never matched by substring, stripped accession version, or assumed
numeric suffix. Candidate FASTA headers are `genome|protein_id`; normalized per-genome
protein FASTAs preserve original protein IDs. GenBank translations and supplied FASTA
sequences remain authoritative; GFF identity alone cannot prove annotation/sequence
concordance, so supply matched files.

### Neighborhood semantics

The default is **eight CDSs on each genomic side**, ordered by `(start,end,ID)`
within the same contig, counting every CDS, including unannotated or untranslated
ones. Order is independent of FASTA or GFF row order. Tied starts are broken
lexicographically and cannot by themselves establish transcriptional relationships.

- `--genes 12`: twelve CDSs on each side.
- `--bp 10000`: also require an inter-CDS gap of at most 10 kb (intersection
  with the gene window).
- `--genes 0 --bp 10000`: use only a base-pair window.
- Coordinates are **1-based inclusive**. Adjacent CDSs have gap 0. Positive
  separation means intervening bases, negative separation means overlap;
  `overlap_bp` gives the intersection length. `intervening_cds = abs(step)-1`.
- `signed_cds_step` is genomic direction, not transcription direction.
  `relative_to_mhc_transcription` accounts for the MHC strand. Orientation is
  co-directional, convergent, divergent, or unknown; none implies an operon.
- All nearby CDSs are included in `neighbors.tsv`, with `qualifies_surface_association`
  distinguishing domain-positive associations from ordinary context. Untranslated
  neighbors are retained with `neighbor_has_protein=False` and cannot have a protein hit.
- Linear boundaries never wrap. Gene-window truncation and unknown/partial anchor
  completeness are flagged. GFF `partial` flags are retained; absence means unknown.
- Circularity is used only from GenBank `topology=circular`, GFF `region` with
  `Is_circular=true`, or an explicit `--circular-contigs` TSV:

  ```tsv
  genome	contig	length
  sampleA	chromosome	2500000
  ```

  For circular contigs the shortest CDS-step path is reported once per neighbor
  (ties favor increasing coordinates). Spacing follows that path; the bp window
  intersects it. This same path definition is used in bp-only mode. Diagrams
  unroll wrapped coordinates and preserve nucleotide spacing. On very small
  circular contigs, eight CDSs on either side may revisit the same genes; output
  never duplicates them.
- Multipart CDSs are mapped and flagged. Neighborhoods on a contig containing
  a compound/split CDS are conservatively **disabled**, because a bounding box
  can span other genes or the replication origin. The protein analysis still
  finishes. This is a documented limitation, not a claimed reconstruction of
  complex gene order. Flatten/curate that annotation only if scientifically justified.

## Tool contributions and evidence interpretation

| Source | What is used | What is not inferred |
|---|---|---|
| MHCscan | Shared original CXXCH+CXXXCH characterization, protein/CDS context concept; new coordinate-safe neighborhood engine | Legacy cluster filtering is not applied to the broad candidate set |
| FeGenie | Existing iron_reduction and iron_oxidation HMMs and per-model full-sequence bit-score cutoffs, searched against all proteins with HMMER | Full FeGenie cluster classification, DIAMOND validation and phenotype calls are not run |
| FEET | Published 3-canonical/5-total motif screen, optional execution of original motif script, imported CELLO/localization, same-contig ±8 CDS concept, nearby MHC reporting | Gram-negative porin/outer-membrane eligibility is not required or used as a gate |
| Pfam | Gathering-threshold domain search on all proteins, including repeated hits | Surface exposure, EET activity, interaction with an adjacent MHC, or sortase substrate identity |

The legacy `bin/MHCscan.sh` contains absolute author-machine paths and runs
Phobius before filtering clusters. Its optional HHsuite annotation step is
commented out. The wrapper therefore reuses the motif helper through a shared
module and imports optional topology data instead of invoking that shell script.
Only `cytoscan.py`'s motif-count expression was refactored; its legacy count is
unchanged. No Phobius binaries are invoked or newly redistributed.

FeGenie's `hmms/iron` directory is preferred (matching its installation setup),
with a fallback to `iron` for other layouts. Only the selected top-level category
HMMs are searched; alignments/reference FASTAs are not mistaken for models.
`--fegenie-categories iron_reduction` restricts the category list. Fractional
cutoffs are preserved; the inspected legacy shell command formats them with
`%d`, truncating fractional values. This wrapper intentionally uses the full
curated numeric cutoff. Every model and threshold is recorded. Both iron
oxidation and reduction hits are retained as family evidence with their source
category, not assertions of reaction direction in the organism.

### Motif definitions and FEET discrepancy

The initial screen is `--min-motifs 2` canonical **CXXCH**, counted at every
starting position. Positions, protein length and motif density are reported.
MHCscan's historical count is non-overlapping **CXXCH + CXXXCH** and is exported
separately as `mhcscan_legacy_motif_count`.

Olmsted et al., *Environmental predictors of electroactive bacterioplankton in
small boreal lakes*, DOI [10.1111/1462-2920.16314](https://doi.org/10.1111/1462-2920.16314),
reports at least **3 CXXCH** and **5 total** motifs, using CXCH, CXXCH, CXXXCH and
a long noncanonical motif. This is a motif eligibility criterion, not the full
published localization/PCC criterion.

The inspected FEET script uses exactly:

```text
C..CH
C.CH
C...CH
C...........[!=C]..CH
```

Python `[!=C]` means a character from `!`, `=`, or `C`; it does **not** mean
"not C". The paper prints the same unusual inequality notation. We do not
silently fix it or assert that the authors intended one particular correction.
Instead the outputs distinguish:

- `feet_code_counts`, `feet_code_total`, `feet_code_motif_pass`: literal,
  non-overlapping `re.findall` compatibility with the source script, independently
  counted per regex. Cross-pattern overlaps are possible; totals are not observed hemes.
- `feet_intended_long_count`, `feet_intended_total`, `feet_intended_motif_pass`:
  an explicitly labeled interpretation replacing `[!=C]` with `[^C]`:
  `C.{11}[^C].{2}CH`. This is a sensitivity comparison, not a claimed author correction.
- `feet_original_script_rows`: unmodified original-script output when `--feet-repo`
  is supplied. This runs its motif script with `3 5`; it does not run the complete
  manual FEET workflow. The external checkout's commit is recorded.
- `feet_full_pipeline_status=not_assessed_Gram_positive_adaptation`: the original
  outer-membrane/periplasm/porin rule is never silently equated with the adaptation.

### Evidence categories, not a confidence score

An MHC may carry `same_protein:binding_domain`, `same_protein:binding_repeat`,
`same_protein:surface_context`, `neighbor:binding_domain`,
`neighbor:anchoring_machinery`, `FeGenie_family_evidence`, and/or a conservative
`LPXTG_C_terminal_architecture_heuristic`. These labels can coexist. Source rows
are available alongside labels. A sortase HMM on the MHC itself is labeled
`same_protein:anchoring_machinery`, not substrate anchoring.

Do not equate a neighboring binding protein with a binding domain in the MHC.
Proximity does not establish an operon, physical interaction, extracellular
exposure, or EET. Conversely, absence of these six domains cannot exclude
other attachment mechanisms. All MHC candidates remain in output even without
qualifying neighbors or surface evidence. `surface_domain_status` and
`fegenie_status` distinguish searches from imports and unassessed analyses.

## Surface databases

The defaults below were verified on **2026-11-20**, via authoritative InterPro
API metadata and the downloadable official Pfam HMM headers. Official Pfam
release metadata reported **39.0**, **2026-09**, UniProtKB **2026_02**. All six
accessions returned current entries and available HMMs; none required substitution.

| Accession/model version | Current short name / name | GA sequence / domain | Interpretation and limitations |
|---|---|---|---|
| PF01473.28 | Choline_bind_1 / Putative cell wall binding repeat | 20.5 / 9.8 | Repeat associated with choline/cell-wall recognition in some proteins; ligands vary, boundaries are approximate |
| PF04122.19 | CW_binding_2 / Cell wall binding domain 2 (CWB2) | 24.9 / 22.3 | Cell-wall binding; repeat architecture matters; PSII binding established in particular clostridial proteins |
| PF00395.27 | SLH / S-layer homology domain | 20.8 / 14.5 | Association with secondary cell-wall polymers; motif alone does not specify ligand affinity or exposure |
| PF04203.20 | Sortase / Sortase domain | 26.6 / 26.6 | Anchoring transpeptidase **machinery**, not a cell-wall binding module on its substrate |
| PF01476.28 | LysM / LysM domain | 20.9 / 11.8 | Peptidoglycan-associated binding; also found in cell-wall degradation enzymes |
| PF07501.19 | G5 / G5 domain | 23 / 23 | Extracellular, adhesion and cell-wall metabolism contexts; not a universal cell-wall anchor |

Primary sources: [PF01473](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01473/),
[PF04122](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF04122/),
[PF00395](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00395/),
[SLH interpretation IPR001119](https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR001119/),
[PF04203](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF04203/),
[PF01476](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01476/),
[PF07501](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF07501/),
[Pfam release metadata](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz).
The HMM download endpoint appends `?annotation=hmm` to each Pfam API URL.

`setup-db` stores the exact downloaded HMM collection and a JSON audit with
retrieval time, URLs, current accession versions, headers, thresholds, metadata
and SHA256 checksums. Keep the collection and audit together for reproducibility.
A future `setup-db` may obtain newer models; reuse the saved collection to reproduce
a run. The analysis manifest includes model versions, hashes and matching setup
provenance. It refuses to adopt an adjacent audit whose checksum mismatches the HMM.

Customize the family list or extract from a versioned local Pfam release:

```bash
mhc-eet setup-db --output custom_db \
  --accessions PF01473,PF04122,PF00395,PF04203,PF01476,PF07501,PF00746

mhc-eet setup-db --output pinned_db \
  --pfam-hmm /path/to/Pfam-A.hmm.gz
```

A deleted/unavailable model causes a visible error; no family replacement is
guessed. `--surface-db your_models.hmm` also supports user-built HMMs directly.
GA is preferred whenever present. Models without GA require an explicit fallback,
for example `--domain-evalue 1e-5`, applied to both sequence and independent domain
E-values; the choice is recorded. This is not claimed to be calibrated for every
user model. No arbitrary minimum coverage is imposed by default, because real
repeats are short; coverage and all individual domains are exported for review.

## Optional localization and existing-result imports

Create `imports/<genome_stem>/` with any supported files below, then pass
`--evidence-dir imports`. Inputs are per genome to prevent collisions between
identically named proteins in different genomes. Unknown IDs or ambiguous mappings
are errors. Raw imported rows are preserved in `imported_source_rows.json`; files
are hashed in the run manifest. Domain imports must use the same protein sequences;
source hashes cannot establish that history automatically.

| File | Format |
|---|---|
| `fegenie.csv` | FeGenie summary with `orf` and original columns such as `category,genome/assembly,orf,HMM,bitscore,bitscore_cutoff`; split multi-genome files before import |
| `feet_heme_count.txt` | Original FEET motif-count TSV with `gene_oid` |
| `mhcscan.csv` | Original MHCscan CSV with `locus`; cluster separator lines are accepted |
| `phobius.short` | Native Phobius `-short` output: sequence ID, TM count, SP, topology |
| `psortb.tsv` | Native PSORTb terse TSV with `SeqID`, `Final_Localization`, optional `Final_Score`; use Gram-positive prediction mode |
| `signalp.tsv` | SignalP 6 `prediction_results.txt`, including `# ID`, `Prediction`, optional `CS Position` |
| `cello.tsv` | Native CELLO "Results" TSV, including `#SeqName` and `#Most-likely-Location` |
| `localization.tsv` | Normalized TSV described below, for any licensed predictor |
| `domains.tsv` | Normalized domain-hit TSV described below |

### CELLO audit and manual route

FEET's `PipelineScript.txt`, steps 3–5, instruct users to concatenate MHCs,
manually upload to CELLO V2.5, click **Results**, save that output, and copy it back.
Its README explicitly discusses the unautomated CELLO step. It is not a local
CLI dependency and is not submitted to any website by this wrapper.

1. Run the core analysis. `genomes/<genome>/proteins.faa` preserves exact per-genome
   IDs; select candidate IDs from `mhcs.tsv` to create a per-genome submission FASTA.
   If using `candidates.faa`, its IDs are prefixed `genome|`; restore exact original
   IDs before import. For reproducibility, record the chosen organism model.
2. If you choose to use CELLO and have access, submit manually and save its native
   Results table as `imports/<genome>/cello.tsv`, keeping IDs unchanged. Alternatively,
   create `localization.tsv` with exact protein IDs and the model/version specified.
3. Rerun with `--evidence-dir imports --resume`. Changed import hashes trigger
   recomputation while identical HMM searches can be reused.

Legacy FEET CELLO results use Gram-negative labels including periplasm and outer
membrane. Native CELLO imports are explicitly marked as organism-model unspecified;
verify their applicability before interpreting Gram-positive candidates. Missing
CELLO data never blocks the core run. An absent localization prediction is
`unavailable`, not "cytoplasmic", "not secreted", or an invented prediction.

### Normalized localization features

Required tab-separated columns: `protein_id`, `source`, `version`, `feature`, `value`.
Extra columns (score, coordinates, organism model, raw prediction) are retained.
Multiple rows per protein and feature are supported, including conflicting sources:

```tsv
protein_id	source	version	feature	value
protein_17	PSORTb_Gram_positive	3.0.3	localization	Cellwall
protein_17	SignalP	6.0	secretion_signal	SP(Sec/SPI)
protein_17	TopologyPredictor	user_recorded_version	transmembrane_topology	predicted_C_terminal_helix
```

Allowed features: `localization`, `secretion_signal`, `transmembrane_topology`,
`lipoprotein_signal`, `membrane_attachment`, `cell_wall_association`, and
`extracellular_exposure`. These are independent annotations. Imports do not
convert secretion into extracellular exposure or convert membrane localization
into cell-wall binding. Use normalized files for DeepTMHMM or other topology
outputs and report their tool version and the model used. Native imports with
no embedded version explicitly label it as unspecified/imported; use normalized
imports when exact version metadata are known.

The built-in LPXTG heuristic requires the motif within the last 80 aa, a downstream
15-aa window containing at least 11 residues from AILMFWVY, and a terminal tail of
3–15 aa containing at least two K/R residues. It reports coordinates and the
architecture test separately from secretion/localization. An isolated C-terminal
LPXTG is recorded with `architecture_pass=false`; an internal motif is not promoted.
This heuristic is neither a TM topology predictor nor a comprehensive search for
all sortase classes or sorting signals. In particular, noncanonical sortase
substrates can be missed. No signal peptide, lipidation or surface exposure is
invented from this heuristic.

### Normalized domain imports

Required TSV columns: `protein_id,source,accession,name,target_from,target_to,domain_bitscore,independent_evalue,model_coverage,threshold`.
Coordinates are 1-based inclusive and checked against the matched protein length.
Add model version, full-sequence score/E-value, HMM coordinates and provenance as
extra columns where available. One row per domain, including repeats. Known family
roles are attached by base accession; unknown configured families are labeled
`user_domain`. `--skip-surface` permits import-only or sequence-only operation and
records the skipped search. Imported hit lists alone do not establish that every
protein was searched: their `surface_domain_status` says so explicitly.

## Outputs

| Output | Contents |
|---|---|
| `mhcs.tsv` | Every MHC passing the broad threshold: motifs, FEET comparisons, source-specific FeGenie and MHCscan results, same-protein domains, independent localization features, neighbor summaries, limitations |
| `neighbors.tsv` | Every MHC–nearby CDS pair: both IDs, coordinates/strands, CDS step/intervening count, gap/overlap, orientation, circular wrap, domains, FeGenie and localization context |
| `domain_hits.tsv` | Pfam/user and FeGenie domain rows across **all** input proteins, including non-MHCs; source/category, model version, thresholds, scores, E-values, HMM/target/envelope coordinates, coverage |
| `genomes.tsv` | Protein/CDS/candidate counts, evidence summaries and input warnings |
| `candidates.faa` | Candidate sequences, namespaced by genome |
| `localization_evidence.tsv` | All imported localization feature rows with provenance |
| `report.html` | Self-contained offline report, separate searchable candidate/neighbor tables, embedded downloadable result files, coordinate-scaled arrow diagrams with hover annotations |
| `manifest.json` | Parameters, input/import/source hashes, repository commits/dirty status, software/database/model versions, thresholds, commands/API calls, skipped analyses and result checksums |
| `genomes/<genome>/` | Normalized proteins, CDS mapping, predicted GFF3 if needed, raw HMMER domtblouts, FeGenie hit table, original FEET motif files/log, imported source rows, search caches |
| `run.log` | Progress and dependency/input diagnostics |

Structured per-protein fields are JSON strings in TSV cells to retain multiple
hits, repeats and conflicting sources. All output text in HTML is escaped. No
remote JavaScript, images, fonts or internet connection is needed to view a report.
Reports can become large with many candidates because data and diagrams are
embedded. Long crowded gene labels may overlap; hover/TSV provides exact values.

## Resume and failure handling

`--resume` reuses a completed analysis only when source/input/database/import
hashes, versions, parameters and output checksums match. Otherwise it recomputes
summaries. Individual HMM searches have content-keyed caches, so an interrupted
run or a changed localization import can reuse matching completed searches.
Cached search results require matching raw domtblout checksums. `--force` recomputes
without accepting old search caches. Partial progress is recorded in
`manifest.inprogress.json`; failed runs never retain a completion marker from a
previous run. Use a new output directory to retain independent analyses.

Optional analyses not requested/available are listed as skipped. Explicitly
requested searches with missing models, incompatible files or a missing optional
Python package produce an actionable error. They are not silently replaced by
empty hits. The wrapper does not download or run predictors without an explicit
user action.

## Validation

See [VALIDATION.md](VALIDATION.md) for exact tests, characterized reference context,
external analyses that were run, and analyses that were not run.
See [REPOSITORY_AUDIT.md](REPOSITORY_AUDIT.md) for integration and license decisions.
