# EVOLVE Hepatitis B analysis

Analysis code for the EVOLVE/Vukuzazi Hepatitis B seroepidemiology paper.

This repository contains the Stata and R workflows used to construct the study cohort, derive survey weights and hepatitis B serostatus outcomes, fit the paper models, and produce tables, figures, and maps. Participant-level source data are confidential and are not included.

## Repository structure

```text
.
|-- config/
|   `-- data_manifest.csv
|-- data/
|   |-- raw/                 # confidential inputs; not tracked
|   `-- derived/             # generated analysis datasets; not tracked
|-- docs/
|   `-- reproducibility_audit.md
|-- outputs/
|   |-- figures/             # generated figures; not tracked
|   |-- logs/                # generated logs; not tracked
|   `-- tables/              # generated tables; not tracked
`-- src/
    |-- R/
    `-- stata/
```


## Software requirements

- Stata 18
- R and the packages checked by `src/R/check_inputs.R`

The R scripts report missing packages but do not install them automatically. This avoids changing the software environment during an analysis run.

## Data configuration

### Option 1: repository-local data

Place approved input files in `data/raw/` and mapping shapefiles in `data/raw/shapefiles/`. The required filenames are listed in `config/data_manifest.csv`.

### Option 2: data stored outside the repository

Set these environment variables before running the analysis:

```sh
export EVOLVE_HBV_ROOT="/path/to/evolve-hepatitis-b"
export EVOLVE_HBV_DATA_RAW="/path/to/approved/input-data"
export EVOLVE_HBV_OUTPUTS_SOURCE="/path/to/existing/model-outputs"
export EVOLVE_HBV_SHAPEFILES="/path/to/shapefiles"
```

`EVOLVE_HBV_DATA_RAW` is used by both Stata and R. `EVOLVE_HBV_OUTPUTS_SOURCE` and `EVOLVE_HBV_SHAPEFILES` are required only by R workflows that consume existing outputs or mapping data.

For an interactive Stata session, the raw-data path may instead be supplied as a global before invoking the runner:

```stata
global EVOLVE_HBV_DATA_RAW "/path/to/approved/input-data"
do "/path/to/evolve-hepatitis-b/src/stata/00_run_all.do"
```

## Running the analysis

### Full Stata pipeline

`src/stata/00_run_all.do` is the canonical Stata entry point. It resolves the repository root from its own filename, runs the input preflight, and then executes the numbered analysis scripts in order.

From Stata:

```stata
do "/path/to/evolve-hepatitis-b/src/stata/00_run_all.do"
```

The runner executes:

1. `08_survey_weights.do` - constructs survey weights.
2. `01_phase1_main_analysis.do` - builds the primary analysis dataset.
3. `02_diagram_numbers.do` - produces participant-flow counts.
4. `03_immunity_threshold_sensitivity.do` - evaluates immunity thresholds.
5. `04_threshold_10.do` - fits the five paper-table survey logistic models.
6. `05_immunity_analysis.do` - prepares immunity categories and plotting data.
7. `06_discordance_analysis.do` - evaluates HBsAg/anti-HBc discordance.
8. `07_treatment_eligibility.do` - derives treatment-eligibility measures.

To check inputs without running the full analysis:

```stata
do "/path/to/evolve-hepatitis-b/src/stata/00_preflight.do"
```

### R workflows

Run these commands from the repository root:

```sh
Rscript src/R/check_inputs.R
Rscript src/R/plot_immunity_figure.R
Rscript src/R/seroepi_graphs.R
Rscript src/R/mapping_evolve.R
```

The input check returns a nonzero status when core files or packages are missing. Mapping dependencies and shapefiles are checked separately.

## Paper model outcomes

`src/stata/04_threshold_10.do` fits survey-weighted logistic regression models for:

1. HBV infection: HBsAg signal-to-cutoff ratio (S/CO) `> 1`.
2. HBV infection: HBsAg S/CO `>= 10`.
3. Vaccine-mediated immunity: HBsAg S/CO `< 10`, anti-HBs positive, and anti-HBc negative.
4. HBV exposure and clearance: HBsAg S/CO `< 10` and anti-HBc positive.
5. Susceptible: HBsAg S/CO `< 10`, anti-HBc negative, and anti-HBs negative.

The internal names `vax_med_immune_ge10`, `hbvexpo_clear_ge10`, and `susceptible_ge10` indicate that an HBsAg threshold of 10 defines infection; consequently, participants assigned to the three non-infection outcomes have HBsAg S/CO `< 10`.

The adjusted models include  HBV vaccine age epoch, sex, education, socioeconomic status, alcohol intake, HIV status, BMI category, hypertension, diabetes, and smoking.

Smoking is modelled as:

- `0`: never smoked (`smoke3cat == 3`)
- `1`: ever smoked, current or former (`smoke3cat == 1` or `smoke3cat == 2`)

## Outputs

Generated files are written to:

- `data/derived/` for intermediate and analysis-ready datasets
- `outputs/tables/` for tables
- `outputs/figures/` for figures and maps
- `outputs/logs/` for Stata logs

These directories are excluded from version control by default. Only disclosure-approved, non-identifiable outputs should be added to a public release.

## Reproducibility safeguards

- Shared Stata and R configuration files centralize paths.
- Preflight scripts stop early when required inputs are unavailable.
- Stata scratch datasets use `tempfile` rather than fixed temporary filenames.
- Derived variables are dropped or safely replaced before recreation, preventing stale variables from contaminating repeated runs.
- Generated Stata variables matching `__0*` are removed before relevant exports.
- Raw data, derived participant-level data, logs, and generated outputs are ignored by Git.

## Data availability

The participant-level datasets cannot be distributed through this repository. Data access is subject to the relevant AHRI governance, ethical approvals, and data-access procedures.

## Citation

Paper citation and DOI will be added after publication. Until then, please cite this repository together with the associated EVOLVE/Vukuzazi Hepatitis B manuscript.

## Contact

For questions about the analysis, open a GitHub issue or contact the study authors through the correspondence details in the paper.

