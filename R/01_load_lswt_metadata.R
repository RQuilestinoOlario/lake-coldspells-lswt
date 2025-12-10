#!/usr/bin/env Rscript
# ================================================================
# LSWT data loading and cold-spell sanity checks
# ------------------------------------------------
# Description: 
#   - Load lake surface water temperature (LSWT) data from HDF5/MAT files
#   - Inspect lake and date metadata
#   - Extract example time series for plotting
#   - Run a small cold-spell (lake cold-spell) test using heatwaveR
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

# Optional: increase memory limit (platform-dependent; may have no effect on some systems)
# mem.maxVSize(320000000000)

# ------------------------------------------------
# 1. Packages
# ------------------------------------------------

## Installation notes (run manually as needed):
## if (!requireNamespace("rhdf5", quietly = TRUE)) {
##   if (!requireNamespace("BiocManager", quietly = TRUE)) {
##     install.packages("BiocManager")
##   }
##   BiocManager::install("rhdf5")
## }
## install.packages(c("dplyr", "tidyr", "ggplot2", "purrr", "heatwaveR"))

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(heatwaveR)
})

# ------------------------------------------------
# 2. File paths
# ------------------------------------------------

# Directory containing the LSWT input files
data_dir <- "your/path/to/data"

# Input file paths
f_lake <- file.path(data_dir, "daily_LSWT_lakeinfo.mat")
f_date <- file.path(data_dir, "daily_LSWT_dateinfo.mat")
f_temp <- file.path(data_dir, "daily_LSWT_data.mat")

# Simple existence checks
stopifnot(
  "Lake info file not found."  = file.exists(f_lake),
  "Date info file not found."  = file.exists(f_date),
  "Temperature file not found" = file.exists(f_temp)
)

# ------------------------------------------------
# 3. Lake metadata
# ------------------------------------------------

# Inspect structure
h5ls(f_lake)

lakeinfo_mat <- h5read(f_lake, "daily_LSWT_lakeinfo")
dim(lakeinfo_mat)

lakeinfo <- as.data.frame(lakeinfo_mat)
colnames(lakeinfo) <- c("Hylak_id", "lat", "lon")

head(lakeinfo)
summary(lakeinfo)

# ------------------------------------------------
# 4. Date metadata
# ------------------------------------------------

h5ls(f_date)

date_mat <- h5read(f_date, "daily_LSWT_dateinfo")
dim(date_mat)

dateinfo <- as.data.frame(date_mat)
colnames(dateinfo) <- c("year", "month", "day")

dateinfo$date <- as.Date(sprintf(
  "%04d-%02d-%02d",
  dateinfo$year,
  dateinfo$month,
  dateinfo$day
))

head(dateinfo)
range(dateinfo$date)

# ------------------------------------------------
# 5. LSWT data: quick peek
# ------------------------------------------------

h5ls(f_temp)

# First 100 lakes, first 365 days
lswt_sub <- h5read(
  file  = f_temp,
  name  = "daily_LSWT_data",
  index = list(1:100, 1:365)
)

dim(lswt_sub)
range(lswt_sub, na.rm = TRUE)

# ------------------------------------------------
# 6. Example time series for one lake
# ------------------------------------------------

lake_i <- 1  # first lake in the file

ts_test <- tibble::tibble(
  date = dateinfo$date[1:365],
  lswt = as.numeric(lswt_sub[lake_i, ])
)

head(ts_test)

ggplot(ts_test, aes(x = date, y = lswt)) +
  geom_line() +
  labs(
    x     = NULL,
    y     = "LSWT (K)",
    title = paste("Example lake", lakeinfo$Hylak_id[lake_i])
  )

# Convert from Kelvin to Celsius
ts_test_c <- tibble::tibble(
  date   = dateinfo$date[1:365],
  lswt_k = as.numeric(lswt_sub[lake_i, ]),
  lswt_c = as.numeric(lswt_sub[lake_i, ]) - 273.15
)

ggplot(ts_test_c, aes(x = date, y = lswt_c)) +
  geom_line() +
  labs(
    x     = NULL,
    y     = "LSWT (°C)",
    title = paste("Example lake", lakeinfo$Hylak_id[lake_i])
  )

# Also keep a Celsius matrix for later use if desired
lswt_sub_c <- lswt_sub - 273.15

# ------------------------------------------------
# 7. Cold-spell test for one lake (heatwaveR)
# ------------------------------------------------

# Settings for the test period / cold-spell detection
test_years         <- 1981:1985
climatology_start  <- "1981-01-01"
climatology_end    <- "1985-12-31"
coldspell_pctile   <- 10
window_half_width  <- 5
min_duration_days  <- 5

idx_test <- which(format(dateinfo$date, "%Y") %in% test_years)

# Load one lake over the selected years and convert to °C
lswt_1 <- h5read(
  file  = f_temp,
  name  = "daily_LSWT_data",
  index = list(1, idx_test)
) - 273.15

ts_lake1 <- tibble::tibble(
  t    = dateinfo$date[idx_test],
  temp = as.numeric(lswt_1)
)

head(ts_lake1)

clim_lake1 <- ts2clm(
  data               = ts_lake1,
  x                  = t,
  y                  = temp,
  climatologyPeriod  = c(climatology_start, climatology_end),
  pctile             = coldspell_pctile,
  windowHalfWidth    = window_half_width
)

lcs_lake1 <- detect_event(
  data        = clim_lake1,
  coldSpells  = TRUE,
  min_duration = min_duration_days
)

str(lcs_lake1$event, max.level = 1)
head(lcs_lake1$event)

# Quick visualization
event_line(
  lcs_lake1,
  spread      = 200,
  metric      = intensity_cumulative,
  start_date  = climatology_start,
  end_date    = climatology_end
)

lolli_plot(
  lcs_lake1,
  metric = intensity_cumulative,
  xaxis  = event_no
)

# ------------------------------------------------
# 8. Cold-spell summary for a few lakes
# ------------------------------------------------

lake_ids_test <- 1:5

lcs_summary_test <- purrr::map_dfr(lake_ids_test, function(i) {
  lswt_i <- h5read(
    file  = f_temp,
    name  = "daily_LSWT_data",
    index = list(i, idx_test)
  ) - 273.15
  
  ts_i <- tibble::tibble(
    t    = dateinfo$date[idx_test],
    temp = as.numeric(lswt_i)
  )
  
  clim_i <- ts2clm(
    data               = ts_i,
    x                  = t,
    y                  = temp,
    climatologyPeriod  = c(climatology_start, climatology_end),
    pctile             = coldspell_pctile,
    windowHalfWidth    = window_half_width
  )
  
  lcs_i <- detect_event(
    data         = clim_i,
    coldSpells   = TRUE,
    min_duration = min_duration_days
  )
  
  # Yearly frequency summary for this lake
  dplyr::count(
    lcs_i$event,
    year = format(date_start, "%Y"),
    name = "n_events"
  ) |>
    dplyr::mutate(lake_idx = i)
})

lcs_summary_test