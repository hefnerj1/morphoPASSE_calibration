[README.md](https://github.com/user-attachments/files/26606782/README.md)
# morphoPASSE calibration

This repository contains the analysis for a study evaluating whether posterior probabilities produced by **MorphoPASSE** are calibrated to the discriminatory capacity of cranial and pelvic morphoscopic traits.

The repository is organized for reproducibility and public release, but the **manuscript is intentionally not included yet**. It will be added after publication.

**Authors: Joseph T. Hefner and Kenzie A. Burns**

## Contents

- a build script to assemble the analysis dataset from the source Excel workbook
- a cleaned analysis script used to generate the current multi-panel Figure 1
- the source workbook used to create `dat_clean`
- the current exported Figure 1 file

## What this repository does not contain

- the manuscript text
- supplemental materials not yet finalized for release


## Repository structure

```text
morphopasse-calibration/
├── README.md
├── LICENSE
├── CITATION.cff
├── code/
│   ├── 00_build_dat_clean.R
│   ├── 01_data_analysis.R
│   ├── README.md
│   └── archive/
├── data/
│   ├── README.md
│   ├── raw/
│   └── derived/
├── figures/
│   ├── README.md
│   ├── generated/
│   └── reference/
└── output/
```

## Analysis workflow

Run the scripts in this order from the repository root.

### 1. Build the derived dataset

`code/00_build_dat_clean.R`

This script:
- reads `data/raw/raw_Data_morpho.xlsx`
- identifies sheet-level dataset and region labels
- combines sheets into a unified analysis object
- writes:
  - `data/derived/dat_clean.rds`
  - `data/derived/dat_clean.csv`

### 2. Run the main analysis

`code/01_data_analysis.R`

This script:
- loads `data/derived/dat_clean.rds`
- identifies morphoscopic trait columns
- computes permutation-based Random Forest benchmarks
- summarizes MorphoPASSE posterior probabilities
- produces the current multi-panel Figure 1
- writes:
  - `figures/generated/Figure1_MultiPanel_A-C.png`

## Getting started
A typical local run looks like this:

```r
source("code/00_build_dat_clean.R")
source("code/01_data_analysis.R")
```

## Citation and license

- Citation metadata are stored in `CITATION.cff`
- Repository reuse is governed by the included `LICENSE`

## Planned additions

After publication:
- the manuscript
- supplemental figures and tables
- additional documentation for peer review and reuse
