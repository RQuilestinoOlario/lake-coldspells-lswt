#!/usr/bin/env Rscript
# ================================================================
# Trend-attribution ratio (TAR) decomposition and maps (Figure 7)
# ------------------------------------------------
# Description:
#   Part A – Detrended LCS metrics (internal variability):
#     - Detrend LSWT time series per lake (linear trend vs time).
#     - Detect cold-spells (LCS) on detrended series (internal variability, IV).
#     - Save annual LCS metrics by lake and year.
#
#   Part B – Observed vs IV vs warming component (ΔSST):
#     - Load observed metrics (from 1981–2010 baseline run) and IV metrics.
#     - Build complete annual series per lake (1981–2020).
#     - Compute warming component as obs − IV for each metric.
#
#   Part C – Trend-attribution ratio (TAR):
#     - For each lake and metric (frequency, duration, intensity),
#       fit linear trends in IV and warming components.
#     - Compute TAR following Oliver-style:
#           TAR = (|T_dSST| − |T_IV|) / max(|T_dSST|, |T_IV|)
#
#   Part D – Figures:
#     - Global time series (Observed / Internal variability / Warming)
#       for frequency, duration, intensity.
#     - Global maps of TAR for frequency, duration, intensity.
#     - Assemble into Figure 7.
#
# Requirements:
#   - LSWT data in HDF5/MAT format: "daily_LSWT_data.mat"
#   - Metadata objects in memory:
#       lakeinfo: Hylak_id, lat, lon, ...
#       dateinfo: date (one row per day; 1981–01–01 to 2020–12–31)
#   - Observed LCS annual metrics already computed and stored in:
#       dir_out = "your/path/to/data/results_lcs_clim1981-2010"
#     (produced by 02_lcs_baseline_sensitivity_figure8.R).
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(heatwaveR)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(maps)
  library(grid)   # for unit()
})

# ------------------------------------------------
# 1. Paths and basic settings
# ------------------------------------------------

# Data directory and LSWT file
data_dir <- "your/path/to/data"
f_temp   <- file.path(data_dir, "daily_LSWT_data.mat")
stopifnot("Temperature file not found." = file.exists(f_temp))

# Index for 1981-01-01 to 2020-12-31
idx_all <- seq_len(nrow(dateinfo))

# Climatology baseline for LCS definition
clim_period <- c("1981-01-01", "2010-12-31")

# Output directory for detrended (IV) metrics
dir_IV  <- file.path(data_dir, "results_lcs_clim1981-2010_detrended")
dir_out <- file.path(data_dir, "results_lcs_clim1981-2010")  # observed (from script 02)

dir.create(dir_IV, showWarnings = FALSE, recursive = TRUE)

# Figure directory
fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Chunking over lakes
n_lakes    <- nrow(lakeinfo)
chunk_size <- 50
lake_ids   <- seq_len(n_lakes)
chunks     <- split(lake_ids, ceiling(seq_along(lake_ids) / chunk_size))

# Base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long)

# ------------------------------------------------
# 2. Detrended LSWT and LCS metrics (internal variability)
# ------------------------------------------------

process_one_lake_IV <- function(lake_idx, date_index = idx_all) {
  
  # 1) Read LSWT (K -> °C)
  lswt_i <- h5read(
    file  = f_temp,
    name  = "daily_LSWT_data",
    index = list(lake_idx, date_index)
  ) - 273.15
  
  ts_i <- tibble::tibble(
    t    = dateinfo$date[date_index],
    temp = as.numeric(lswt_i)
  )
  
  # 2) Linear detrending (internal variability)
  t_num <- as.numeric(ts_i$t)
  fit   <- lm(temp ~ t_num, data = ts_i)
  
  ts_i$temp_IV <- resid(fit)  # detrended series
  
  # 3) Climatology + LCS detection on detrended series
  clim_i <- ts2clm(
    data              = ts_i,
    x                 = t,
    y                 = temp_IV,
    climatologyPeriod = clim_period,
    pctile            = 10,
    windowHalfWidth   = 5
  )
  
  lcs_i <- detect_event(
    data         = clim_i,
    coldSpells   = TRUE,
    min_duration = 5
  )
  
  ev <- lcs_i$event
  
  # 4) No events -> empty tibble
  if (nrow(ev) == 0) {
    return(tibble::tibble(
      year                = integer(),
      n_events            = integer(),
      duration_total      = numeric(),
      duration_mean       = numeric(),
      intensity_mean_mean = numeric(),
      intensity_cum_total = numeric(),
      lake_idx            = integer(),
      Hylak_id            = integer()
    ))
  }
  
  # Ensure Date class
  if (!inherits(ev$date_start, "Date")) {
    ev$date_start <- as.Date(ev$date_start, origin = "1970-01-01")
  }
  
  ev %>%
    dplyr::mutate(year = as.integer(format(date_start, "%Y"))) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      n_events            = dplyr::n(),
      duration_total      = sum(duration,            na.rm = TRUE),
      duration_mean       = mean(duration,           na.rm = TRUE),
      intensity_mean_mean = mean(intensity_mean,     na.rm = TRUE),
      intensity_cum_total = sum(intensity_cumulative, na.rm = TRUE),
      .groups             = "drop"
    ) %>%
    dplyr::mutate(
      lake_idx = lake_idx,
      Hylak_id = lakeinfo$Hylak_id[lake_idx]
    )
}

# Loop over chunks
for (b in seq_along(chunks)) {
  block_ids <- chunks[[b]]
  
  message(
    "IV block ", b, " / ", length(chunks),
    " (lakes ", min(block_ids), "-", max(block_ids), ")"
  )
  
  block_res <- purrr::map_dfr(block_ids, process_one_lake_IV)
  
  out_file <- file.path(
    dir_IV,
    sprintf("lcs_annual_metrics_IV_block_%03d.rds", b)
  )
  saveRDS(block_res, out_file)
}

# ------------------------------------------------------------
# 3. Load annual LCS metrics (observed + detrended = IV)
# ------------------------------------------------------------

files_obs <- list.files(dir_out, pattern = "\\.rds$", full.names = TRUE)
lcs_obs   <- purrr::map_dfr(files_obs, readRDS)

files_IV  <- list.files(dir_IV,  pattern = "\\.rds$", full.names = TRUE)
lcs_IV    <- purrr::map_dfr(files_IV,  readRDS)

# Attach lat/lon via Hylak_id
lake_xy <- lakeinfo %>%
  dplyr::select(Hylak_id, lat, lon)

lcs_obs <- lcs_obs %>%
  dplyr::left_join(lake_xy, by = "Hylak_id") %>%
  dplyr::rename(
    freq_obs = n_events,
    dur_obs  = duration_total,
    int_obs  = intensity_mean_mean
  )

lcs_IV <- lcs_IV %>%
  dplyr::left_join(lake_xy, by = "Hylak_id") %>%
  dplyr::rename(
    freq_IV = n_events,
    dur_IV  = duration_total,
    int_IV  = intensity_mean_mean
  )

# ------------------------------------------------------------
# 4. Complete annual series per lake and ΔSST components
# ------------------------------------------------------------

years_full  <- 1981:2020
ids_tracked <- sort(unique(c(lcs_obs$Hylak_id, lcs_IV$Hylak_id)))

skeleton <- tidyr::expand_grid(
  Hylak_id = ids_tracked,
  year     = years_full
)

annual_pair <- skeleton %>%
  # Observed metrics
  dplyr::left_join(
    lcs_obs %>%
      dplyr::select(
        Hylak_id, year,
        lake_idx_obs = lake_idx,
        lat_obs      = lat,
        lon_obs      = lon,
        freq_obs, dur_obs, int_obs
      ),
    by = c("Hylak_id", "year")
  ) %>%
  # Internal variability metrics
  dplyr::left_join(
    lcs_IV %>%
      dplyr::select(
        Hylak_id, year,
        lake_idx_IV = lake_idx,
        lat_IV      = lat,
        lon_IV      = lon,
        freq_IV, dur_IV, int_IV
      ),
    by = c("Hylak_id", "year")
  ) %>%
  # Coalesce coordinates/indices
  dplyr::mutate(
    lake_idx = dplyr::coalesce(lake_idx_obs, lake_idx_IV),
    lat      = dplyr::coalesce(lat_obs,      lat_IV),
    lon      = dplyr::coalesce(lon_obs,      lon_IV)
  ) %>%
  dplyr::select(
    Hylak_id, lake_idx, lat, lon, year,
    freq_obs, dur_obs, int_obs,
    freq_IV,  dur_IV,  int_IV
  ) %>%
  # Fill missing years with zeros for metrics
  dplyr::mutate(
    freq_obs = tidyr::replace_na(freq_obs, 0),
    dur_obs  = tidyr::replace_na(dur_obs,  0),
    int_obs  = tidyr::replace_na(int_obs,  0),
    freq_IV  = tidyr::replace_na(freq_IV,  0),
    dur_IV   = tidyr::replace_na(dur_IV,   0),
    int_IV   = tidyr::replace_na(int_IV,   0),
    # ΔSST component = obs − IV
    freq_dSST = freq_obs - freq_IV,
    dur_dSST  = dur_obs  - dur_IV,
    int_dSST  = int_obs  - int_IV
  )

# ------------------------------------------------------------
# 5. Helper to compute TAR for a single metric
# ------------------------------------------------------------

compute_tar_one_metric <- function(df, metric_IV, metric_dSST) {
  
  # Require enough time variation
  if (dplyr::n_distinct(df$year) < 10) {
    return(tibble::tibble(
      T_IV   = NA_real_,
      T_dSST = NA_real_,
      TAR    = NA_real_
    ))
  }
  
  lm_IV   <- lm(df[[metric_IV]]   ~ df$year)
  lm_dSST <- lm(df[[metric_dSST]] ~ df$year)
  
  T_IV   <- coef(lm_IV)[2]
  T_dSST <- coef(lm_dSST)[2]
  
  # If both trends essentially zero, TAR undefined
  if (all(abs(c(T_IV, T_dSST)) < 1e-8)) {
    return(tibble::tibble(
      T_IV   = NA_real_,
      T_dSST = NA_real_,
      TAR    = NA_real_
    ))
  }
  
  TAR <- (abs(T_dSST) - abs(T_IV)) / max(abs(T_dSST), abs(T_IV))
  
  tibble::tibble(T_IV = T_IV, T_dSST = T_dSST, TAR = TAR)
}

# ------------------------------------------------------------
# 6. TAR per lake for frequency, duration, intensity
# ------------------------------------------------------------

tar_freq <- annual_pair %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::group_modify(~ compute_tar_one_metric(
    .x,
    metric_IV   = "freq_IV",
    metric_dSST = "freq_dSST"
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    T_IV_decade   = T_IV   * 10,
    T_dSST_decade = T_dSST * 10
  )

tar_dur <- annual_pair %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::group_modify(~ compute_tar_one_metric(
    .x,
    metric_IV   = "dur_IV",
    metric_dSST = "dur_dSST"
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    T_IV_decade   = T_IV   * 10,
    T_dSST_decade = T_dSST * 10
  )

tar_int <- annual_pair %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::group_modify(~ compute_tar_one_metric(
    .x,
    metric_IV   = "int_IV",
    metric_dSST = "int_dSST"
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    T_IV_decade   = T_IV   * 10,
    T_dSST_decade = T_dSST * 10
  )

# ------------------------------------------------------------
# 7. Maps of TAR for freq / duration / intensity
# ------------------------------------------------------------

cols_rdbu <- hcl.colors(11, "RdBu")
cols_rdbu <- cols_rdbu[-c(5, 6)]

fig_TAR_freq <- ggplot(tar_freq %>% dplyr::filter(!is.na(TAR)),
                       aes(x = lon, y = lat)) +
  geom_polygon(
    data         = map_base,
    aes(x = lon, y = lat, group = group),
    inherit.aes  = FALSE,
    colour       = NA,
    fill         = "grey90"
  ) +
  geom_point(aes(colour = TAR), size = 0.4) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lakeinfo$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lakeinfo$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_rdbu),
    limits  = c(-1, 1),
    name    = expression(bold(TAR['freq']))
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.1, 0.4),
    legend.key.height  = unit(0.5, "cm"),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.background  = element_rect(fill = "transparent", colour = NA)
  )

fig_TAR_dur <- ggplot(tar_dur %>% dplyr::filter(!is.na(TAR)),
                      aes(x = lon, y = lat)) +
  geom_polygon(
    data         = map_base,
    aes(x = lon, y = lat, group = group),
    inherit.aes  = FALSE,
    colour       = NA,
    fill         = "grey90"
  ) +
  geom_point(aes(colour = TAR), size = 0.4) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lakeinfo$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lakeinfo$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_rdbu),
    limits  = c(-1, 1),
    name    = expression(bold(TAR['dur']))
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.1, 0.4),
    legend.key.height  = unit(0.5, "cm"),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.background  = element_rect(fill = "transparent", colour = NA)
  )

fig_TAR_int <- ggplot(tar_int %>% dplyr::filter(!is.na(TAR)),
                      aes(x = lon, y = lat)) +
  geom_polygon(
    data         = map_base,
    aes(x = lon, y = lat, group = group),
    inherit.aes  = FALSE,
    colour       = NA,
    fill         = "grey90"
  ) +
  geom_point(aes(colour = TAR), size = 0.4) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lakeinfo$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lakeinfo$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_rdbu),
    limits  = c(-1, 1),
    name    = expression(bold(TAR['int']))
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.1, 0.4),
    legend.key.height  = unit(0.5, "cm"),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.background  = element_rect(fill = "transparent", colour = NA)
  )

# ------------------------------------------------------------
# 8. Global annual time series: Observed vs IV vs Warming
# ------------------------------------------------------------

cols_src <- c(
  "Observed"             = "#E69F00",
  "Internal variability" = "#006DA9",
  "Warming"              = "#AB2C00"
)

pd <- position_dodge(width = 0.7)

# Frequency
global_freq <- annual_pair %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    `Observed`             = mean(freq_obs,  na.rm = TRUE),
    `Internal variability` = mean(freq_IV,   na.rm = TRUE),
    `Warming`              = mean(freq_dSST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = -year,
    names_to  = "source",
    values_to = "value"
  )

fig7a <- ggplot(global_freq, aes(x = year, y = value, fill = source)) +
  geom_col(position = pd, width = 0.55, colour = NA) +
  geom_smooth(
    aes(colour = source, group = source),
    method       = "lm",
    se           = FALSE,
    linewidth    = 0.5,
    inherit.aes  = TRUE,
    show.legend  = FALSE
  ) +
  scale_fill_manual(values = cols_src, name = NULL) +
  scale_colour_manual(values = cols_src, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Year",
    y = "Mean annual LCS\nfrequency (events/year)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.position  = c(0.7, 0.9),
    legend.text      = element_text(size = 10)
  )

# Duration
global_dur <- annual_pair %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    `Observed`             = mean(dur_obs,  na.rm = TRUE),
    `Internal variability` = mean(dur_IV,   na.rm = TRUE),
    `Warming`              = mean(dur_dSST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = -year,
    names_to  = "source",
    values_to = "value"
  )

fig7b <- ggplot(global_dur, aes(x = year, y = value, fill = source)) +
  geom_col(position = pd, width = 0.55, colour = NA) +
  geom_smooth(
    aes(colour = source, group = source),
    method       = "lm",
    se           = FALSE,
    linewidth    = 0.5,
    inherit.aes  = TRUE,
    show.legend  = FALSE
  ) +
  scale_fill_manual(values = cols_src, name = NULL) +
  scale_colour_manual(values = cols_src, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Year",
    y = "Mean annual LCS\nduration (days/year)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# Intensity
global_int <- annual_pair %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    `Observed`             = mean(int_obs,  na.rm = TRUE),
    `Internal variability` = mean(int_IV,   na.rm = TRUE),
    `Warming`              = mean(int_dSST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = -year,
    names_to  = "source",
    values_to = "value"
  )

fig7c <- ggplot(global_int, aes(x = year, y = value, fill = source)) +
  geom_col(position = pd, width = 0.55, colour = NA) +
  geom_smooth(
    aes(colour = source, group = source),
    method       = "lm",
    se           = FALSE,
    linewidth    = 0.5,
    inherit.aes  = TRUE,
    show.legend  = FALSE
  ) +
  scale_fill_manual(values = cols_src, name = NULL) +
  scale_colour_manual(values = cols_src, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Year",
    y = "Mean LCS intensity\n(°C/event)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# ------------------------------------------------------------
# 9. Assemble Figure 7 and save
# ------------------------------------------------------------

fig7 <- plot_grid(
  fig7a, fig_TAR_freq,
  fig7b, fig_TAR_dur,
  fig7c, fig_TAR_int,
  ncol   = 2,
  align  = "v",
  labels = c("(a)", "(d)", "(b)", "(e)", "(c)", "(f)")
)

save_plot(
  filename    = file.path(fig_dir, "fig7_tar.pdf"),
  plot        = fig7,
  base_width  = 14,
  base_height = 8
)