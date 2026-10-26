
# sparse4-hnscc

Reproducibility code for the manuscript *B-cell-associated transcriptional
context and post hoc prognostic modelling in head and neck squamous cell
carcinoma*.

Author and correspondence details are omitted for anonymous peer review.

This repository contains analysis entry points, the locked sparse4 scoring
specification, aggregate result tables, current panel assets that meet the
release size limit, software manifests, random seeds, and a manuscript-to-code
map. It does not redistribute raw data or participant-level derivatives.

## Structure

- `config/`: paths, parameters, locked score material, seeds, and release state.
- `data/`: public accession metadata and local placement instructions.
- `environment/`: R lockfile, Python requirements, and observed versions.
- `R/`: shared validation, I/O, and logging utilities.
- `scripts/`: numbered analysis, scoring, figure, and table entry points.
- `archive_exploratory/`: retained exploratory or secondary historical code.
- `results/`: non-identifying aggregate tables and distributable panel assets.
- `docs/`: scoring specification, limitations, release checklist, and complete figure/table map.
- `provenance/`: source hashes and release checksums.
- `tests/`: static, configuration, and scoring-parity checks.

## Environment

Target versions are R 4.4.1 and Python 3.10.18. Restore R packages with the
lockfile and Python packages with the requirements file. A full rerun requires
Linux, substantial memory for the single-cell/spatial objects, and the public
datasets listed in `data/accessions.csv`. Historical package-version gaps are
described in `docs/known_limitations.md`.

```bash
Rscript -e 'renv::restore(lockfile = "environment/renv.lock")'
python3 -m venv .venv
. .venv/bin/activate
pip install -r environment/requirements.txt
```

## Data acquisition

Follow `data/README.md`. Public accessions include TCGA-HNSC, GSE164690,
GSE208253, GSE41613, GSE65858, and GSE159067. Download dates that were not
recoverable are marked `not recorded`; they were not inferred.

## Running the repository

Start with the fail-closed release checks:

```bash
bash run_all.sh --check-only
```

Test the portable score implementations on a synthetic matrix:

```bash
bash run_all.sh --score-fixture
```

The full analysis dispatcher is available for local environments after the
required public datasets and reconstructed inputs described in `data/README.md`
have been placed under the documented paths:

```bash
bash run_all.sh --full
```

The dispatcher validates the required local inputs and their recorded hashes
before beginning analysis. Historical and secondary procedures are kept
separate under `archive_exploratory/` and are not launched by default.

## Paper correspondence

`docs/figure_table_map.md` and its CSV companion map every current figure panel,
released table, and key numerical claim to the relevant script and evidence.
`results/manifest.csv` records the producing script, paper item, sample count,
seed rule, and creation date for released outputs.
`docs/table_descriptions.md` gives a concise interpretation and evidence
boundary for each main and supplementary table.
The mapped figure renderers preserve provenance, but several still require
non-distributed reconstructed inputs. `scripts/23_make_figures.R` verifies the
released panel assets and map; it is not a universal panel-regeneration engine.

## Randomness and cross-platform variation

Seeds and seed rules are centralized in `config/params.yml` and
`config/random_seeds.tsv`. Parallel random-number use is declared by the
relevant script. Floating-point and solver differences across platforms may
produce negligible numerical variation; they must not change the locked model
identity, gene order, coefficient signs, analysis roster, or evidence-level
interpretation.

## Limitations

Read `docs/known_limitations.md` before interpreting results. In particular,
sparse4 is post hoc, the external cohorts had been inspected during model
development, and the spatial Yellow/TLS analysis is not a spatial validation
of sparse4.

## Licence and citation

Code is released under the MIT License. Documentation is released under
CC BY 4.0. Citation metadata are provided in `CITATION.cff`, and the repository
URL is recorded in `config/release.yml`.
