# Lake cold-spells (LCS) from global lake surface water temperatures

This repository contains R scripts to compute and analyze lake cold-spells (LCS)
from daily lake surface water temperature (LSWT) data, including:

- LSWT data loading and pre-processing
- LCS detection using `heatwaveR`
- Spatial patterns and trends (Figures 1–4)
- Latitudinal and thermal structure (Figure 5)
- Temporal evolution by latitude band (Figure 6)
- Trend-attribution ratio (TAR) decomposition (Figure 7)

## Repository structure

- `R/` – R scripts for data processing and figure generation.
- `LICENSE` – License for this repository.
- `CITATION.cff` – Citation metadata for this repository.

## Data

LSWT and lake metadata are **not included** in this repository due to size and/or licensing.
To reproduce the results, obtain:

- `daily_LSWT_data.mat`
- `daily_LSWT_lakeinfo.mat`
- `daily_LSWT_dateinfo.mat`

Place them in `your/path/to/data` and/or update the `data_dir` paths
in the scripts under `R/`.

## Usage

Example (from R):

```r
# Set working directory to repository root
setwd("your/local/path/lake-coldspells-lswt")

# Run scripts in order:
source("R/01_load_lswt_metadata.R")
source("R/03_lswt_mean_hemispheres_figure1.R")
# ...
