# Repository and dependency audit

Inspected on 2026-11-20 before implementing the extension. Repositories were
resolved through the connected Middle-Author-Bioinformatics GitHub account.

| Repository | Inspected commit | Instructions/source/output inspected |
|---|---|---|
| [Middle-Author-Bioinformatics/FeGenie](https://github.com/Middle-Author-Bioinformatics/FeGenie) | `034f554fb12528b43f7515b5c066de610a6ec6e8` | README, LICENSE, setup scripts, FeGenie.py, HMM library and cutoffs, CSV headers, Thermincola reference FASTA |
| [Middle-Author-Bioinformatics/MHCscan](https://github.com/Middle-Author-Bioinformatics/MHCscan) | `2fe42012d46f3c57281b856331dc51e3777b1872` | README, GPL LICENSE, setup.sh, bin/MHCscan.sh, cytoscan.py, mhc_annotate.py, Phobius script headers |
| [Middle-Author-Bioinformatics/FEET](https://github.com/Middle-Author-Bioinformatics/FEET) | `bd256f4df08fc70b69ede210d75f46466e5d002c` | README, environment/setup/FeGenie-modification instructions, PipelineScript.txt, both motif/neighborhood scripts, example CELLO and pipeline outputs; no LICENSE found |

FEET is a fork of [McMahonLab/FEET](https://github.com/McMahonLab/FEET).
The published methodology was checked against
[Olmsted et al., DOI 10.1111/1462-2920.16314](https://doi.org/10.1111/1462-2920.16314).
The 3 canonical/5 total requirement matches the script parameters. The long-regex
notation discrepancy is documented in the CLI guide and explicitly tested.

## Placement

MHCscan is the host because the new tool starts with a broad MHC screen and extends
protein characterization with genomic context. Its explicit GPLv3 license supports
an extension with a shared motif helper. FeGenie retains its separate database and
source-specific family evidence. FEET's pipeline is partly manual, Gram-negative
oriented, and has no explicit repository license; copying its workflow wholesale
would also carry unsuitable filtering and fragile Prodigal-header assumptions.
The extension runs an explicitly supplied FEET motif script as an independent
program and supports source-preserving result imports. No FEET script/model is
copied into this repository. No original repository is merged automatically.

## Existing functionality and constraints

- **MHCscan:** `cytoscan.py` counts non-overlapping CXXCH and CXXXCH, combines
  Phobius signal information with heme-bearing genes and emits filtered clusters.
  GFF row order or numeric Prodigal suffixes determine its old neighborhoods.
  Its shell wrapper contains absolute machine paths. This extension factors its
  motif-count function into a shared module, retains the legacy count and adds
  a separate canonical threshold and validated coordinate ordering. The legacy
  cluster eligibility filter and optional HHsuite step are not represented as
  newly executed analyses. Legacy CSV imports preserve `locus` and all columns.
- **FeGenie:** uses HMMER model-specific full-sequence cutoffs, Prodigal, clustering,
  optional DIAMOND/BLAST validation, and optional R reporting. This extension uses
  the actual curated HMMs/cutoff files, searches all proteins through PyHMMER,
  and preserves raw domain output and source categories. Its family-search
  integration is not described as full FeGenie pathway classification. Fractional
  cutoff precision is deliberately retained. FeGenie CSV imports require the
  original `orf` field and preserve original source rows.
- **FEET:** its motif script requires BioPython and produces count/list/FASTA
  files. `MHC-BB-OMP-matcher.py` depends on specific Prodigal header fields,
  substring regex matches, numeric ORF indices and an inclusive ±N window.
  The wrapper implements exact-ID, coordinate-ordered same-contig windows with
  all CDSs counted. It reports nearby MHC and surface components without
  requiring a porin. Full FEET uses FeGenie, HMMER, DIAMOND, Prodigal, BLAST, R,
  shell utilities and a manual CELLO submission. Only the motif executable is
  automatically invoked here; CELLO and pre-existing results can be imported.

## License and installation choices

| Component | Inspected license / availability | Decision |
|---|---|---|
| MHCscan | GPL version 3 LICENSE | Extension distributed under GPL-3.0-only alongside existing license |
| FeGenie | AGPL version 3 LICENSE | External checkout/database only, not copied or relicensed; retain original notices and terms |
| FEET | No explicit LICENSE found at inspected commit | No redistribution of its scripts/models; optional external executable or imported results |
| BioPython | BioPython License Agreement | Core parser dependency, install from PyPI |
| PyHMMER | MIT bindings; HMMER/Easel BSD licenses | Optional analysis extra; avoids system HMMER installation |
| Pyrodigal / Prodigal | GPLv3 | Optional analysis extra; deterministic gene prediction without separate executable |
| Pfam / InterPro data | EMBL-EBI data and database terms, attribution required | Download explicitly from official endpoints, preserve provenance; no models committed |
| PSORTb | GPL standalone; Gram-positive mode; Perl/BioPerl and other installation dependencies | Import terse output or normalized features; no automatic heavyweight installation |
| SignalP 6 | Academic standalone download; commercial users require DTU licensing | Import predictions only; no automated download/installation or assumed commercial permission |
| Phobius | Existing repository contains binaries/models; separate distribution terms were not established for a new distribution | Import native short output only; do not invoke or newly bundle it |
| CELLO V2.5 | FEET documents a manual web workflow, not an installable local dependency | Optional manual import; no automated website submission |
| Other topology/localization tools | Tool/model-specific terms vary | Normalized feature import, explicit source/version; no implicit installation |

The wrapper does not assume an academic-only SignalP entitlement. No paid or
restricted predictor is required.

Primary dependency/license references:

- [PyHMMER copyright](https://pyhmmer.readthedocs.io/en/latest/guide/copyright.html)
- [Pyrodigal copyright](https://pyrodigal.readthedocs.io/en/stable/guide/copyright.html)
- [PSORTb installation and GPL](https://psort.org/downloads/INSTALL.html)
- [PSORTb Gram-positive support](https://psort.org/psortb/)
- [SignalP 6 distribution](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
- [SignalP 6 publication/distribution description](https://pmc.ncbi.nlm.nih.gov/articles/PMC9287161/)
- [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/)

Exact runtime package versions, repository commits, dirty state and source hashes
are included in each analysis manifest. A dirty working checkout is explicitly
reported, not represented as the unmodified upstream commit.
