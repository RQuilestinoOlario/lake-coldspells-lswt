#!/usr/bin/env Rscript
# ================================================================
# Lake cold-spell baseline sensitivity and Figure 8
# ------------------------------------------------
# Description:
#   - Detect lake cold-spells (LCS) using two climatology baselines:
#       * 1981–2010
#       * 1991–2020
#   - Aggregate annual LCS metrics per lake
#   - Compare baselines and generate Figure 8
#
# Requirements:
#   - LSWT data in HDF5/MAT format (daily_LSWT_data.mat)
#   - Lake and date metadata: objects `lakeinfo` and `dateinfo`
#     available in the environment, or loaded beforehand, e.g.:
#       lakeinfo <- readRDS("your/path/to/lakeinfo.rds")
#       dateinfo <- readRDS("your/path/to/dateinfo.rds")
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
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(maps)
})

# ------------------------------------------------
# 1. Paths and basic settings
# ------------------------------------------------

# Directory containing the LSWT input files
data_dir <- "your/path/to/data"

# LSWT data file (daily LSWT, K)
f_temp <- file.path(data_dir, "daily_LSWT_data.mat")

stopifnot("Temperature file not found." = file.exists(f_temp))

# Full time index: 1981-01-01 to 2020-12-31
idx_all <- seq_len(nrow(dateinfo))

# Output directories for chunked results
dir_out_8110 <- file.path(data_dir, "results_lcs_clim1981-2010")
dir_out_9120 <- file.path(data_dir, "results_lcs_clim1991-2020")

dir.create(dir_out_8110, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_out_9120, showWarnings = FALSE, recursive = TRUE)

# Figure output directory
fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Chunking over lakes
n_lakes    <- nrow(lakeinfo)
chunk_size <- 50
lake_ids   <- seq_len(n_lakes)

chunks <- split(lake_ids, ceiling(seq_along(lake_ids) / chunk_size))
length(chunks)

# ------------------------------------------------
# 2. Helper: process one lake for a given climatology
# ------------------------------------------------

process_one_lake <- function(lake_idx,
                             date_index   = idx_all,
                             clim_period  = c("1981-01-01", "2010-12-31"),
                             f_temp       = f_temp,
                             dateinfo_obj = dateinfo,
                             lakeinfo_obj = lakeinfo) {
  
  # 1. Read LSWT (K -> °C)
  lswt_i <- h5read(
    file  = f_temp,
    name  = "daily_LSWT_data",
    index = list(lake_idx, date_index)
  ) - 273.15
  
  ts_i <- tibble::tibble(
    t    = dateinfo_obj$date[date_index],
    temp = as.numeric(lswt_i)
  )
  
  # 2. Climatology + LCS detection
  clim_i <- ts2clm(
    data              = ts_i,
    x                 = t,
    y                 = temp,
    climatologyPeriod = clim_period,
    pctile            = 10,
    windowHalfWidth   = 5
  )
  
  lcs_i <- detect_event(
    data        = clim_i,
    coldSpells  = TRUE,
    min_duration = 5
  )
  
  ev <- lcs_i$event
  
  # 3. No events -> empty tibble with correct columns
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
  
  # 4. Ensure Date class for date_start
  if (!inherits(ev$date_start, "Date")) {
    ev$date_start <- as.Date(ev$date_start, origin = "1970-01-01")
  }
  
  # 5. Add year and summarise by year
  ev %>%
    dplyr::mutate(year = as.integer(format(.data$date_start, "%Y"))) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      n_events            = dplyr::n(),
      duration_total      = sum(duration, na.rm = TRUE),
      duration_mean       = mean(duration, na.rm = TRUE),
      intensity_mean_mean = mean(intensity_mean, na.rm = TRUE),
      intensity_cum_total = sum(intensity_cumulative, na.rm = TRUE),
      .groups             = "drop"
    ) %>%
    dplyr::mutate(
      lake_idx = lake_idx,
      Hylak_id = lakeinfo_obj$Hylak_id[lake_idx]
    )
}

# ------------------------------------------------
# 3. Run LCS pipeline for two climatology baselines
# ------------------------------------------------

run_lcs_for_climatology <- function(clim_period,
                                    output_dir,
                                    chunks_list = chunks) {
  message("Running LCS detection for climatology: ",
          clim_period[1], " to ", clim_period[2])
  
  for (b in seq_along(chunks_list)) {
    block_ids <- chunks_list[[b]]
    
    message(
      "Block ", b, " / ", length(chunks_list),
      " (lakes ", min(block_ids), "-", max(block_ids), ")"
    )
    
    block_res <- purrr::map_dfr(
      block_ids,
      process_one_lake,
      clim_period = clim_period
    )
    
    out_file <- file.path(
      output_dir,
      sprintf("lcs_annual_metrics_block_%03d.rds", b)
    )
    saveRDS(block_res, out_file)
  }
}

# 3.1. Baseline 1981–2010
clim_period_8110 <- c("1981-01-01", "2010-12-31")
run_lcs_for_climatology(clim_period_8110, dir_out_8110)

# 3.2. Baseline 1991–2020
clim_period_9120 <- c("1991-01-01", "2020-12-31")
run_lcs_for_climatology(clim_period_9120, dir_out_9120)

# ------------------------------------------------
# 4. Load results and attach geometry
# ------------------------------------------------

# 4.1. Load 1981–2010 baseline
files_8110 <- list.files(
  dir_out_8110,
  pattern    = "^lcs_annual_metrics_block_\\d+\\.rds$",
  full.names = TRUE
)

lcs_8110 <- files_8110 |>
  lapply(readRDS) |>
  dplyr::bind_rows()

dim(lcs_8110)
head(lcs_8110)

# 4.2. Load 1991–2020 baseline
files_9120 <- list.files(
  dir_out_9120,
  pattern    = "^lcs_annual_metrics_block_\\d+\\.rds$",
  full.names = TRUE
)

lcs_9120 <- files_9120 |>
  lapply(readRDS) |>
  dplyr::bind_rows()

dim(lcs_9120)
head(lcs_9120)

# 4.3. Attach lon/lat
lake_index_df <- lakeinfo %>%
  dplyr::mutate(lake_idx = dplyr::row_number())

lcs_8110_geo <- lcs_8110 %>%
  dplyr::left_join(lake_index_df, by = "lake_idx")

lcs_9120_geo <- lcs_9120 %>%
  dplyr::left_join(lake_index_df, by = "lake_idx")

lcs_8110_geo <- lcs_8110_geo %>%
  dplyr::mutate(Hylak_id = dplyr::coalesce(Hylak_id.x, Hylak_id.y)) %>%
  dplyr::select(-Hylak_id.x, -Hylak_id.y)

lcs_9120_geo <- lcs_9120_geo %>%
  dplyr::mutate(Hylak_id = dplyr::coalesce(Hylak_id.x, Hylak_id.y)) %>%
  dplyr::select(-Hylak_id.x, -Hylak_id.y)

# ------------------------------------------------
# 5. Compare baselines: annual means over 1991–2020
# ------------------------------------------------

lcs_8110_mean_9120 <- lcs_8110_geo %>%
  dplyr::filter(year >= 1991, year <= 2020) %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::summarise(
    freq_mean_8110     = mean(n_events, na.rm = TRUE),
    dur_total_8110     = mean(duration_total, na.rm = TRUE),
    int_cum_mean_8110  = mean(intensity_cum_total, na.rm = TRUE),
    .groups            = "drop"
  )

lcs_9120_mean_9120 <- lcs_9120_geo %>%
  dplyr::filter(year >= 1991, year <= 2020) %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::summarise(
    freq_mean_9120     = mean(n_events, na.rm = TRUE),
    dur_total_9120     = mean(duration_total, na.rm = TRUE),
    int_cum_mean_9120  = mean(intensity_cum_total, na.rm = TRUE),
    .groups            = "drop"
  )

baseline_comp <- lcs_8110_mean_9120 %>%
  dplyr::inner_join(
    lcs_9120_mean_9120,
    by = c("Hylak_id", "lake_idx", "lat", "lon")
  ) %>%
  dplyr::mutate(
    freq_diff    = freq_mean_9120    - freq_mean_8110,
    dur_diff     = dur_total_9120    - dur_total_8110,
    int_cum_diff = int_cum_mean_9120 - int_cum_mean_8110
  )

head(baseline_comp)

# ------------------------------------------------
# 6. Figure 8: base map and scatter panels
# ------------------------------------------------

map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long)

# 6.1. Scatter comparisons (b, e, h)
fig8b <- ggplot(baseline_comp, aes(x = freq_mean_8110, y = freq_mean_9120)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(shape = 1, colour = "grey50", alpha = 0.3, size = 1) +
  labs(
    x = "Mean annual LCS frequency\n(baseline 1981–2010)",
    y = "Mean annual LCS frequency\n(baseline 1991–2020)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

fig8e <- ggplot(baseline_comp, aes(x = dur_total_8110, y = dur_total_9120)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(shape = 1, colour = "grey50", alpha = 0.3, size = 1) +
  labs(
    x = "Mean annual LCS duration\n(baseline 1981–2010)",
    y = "Mean annual LCS duration\n(baseline 1991–2020)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

fig8h <- ggplot(baseline_comp, aes(x = int_cum_mean_8110, y = int_cum_mean_9120)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(shape = 1, colour = "grey50", alpha = 0.3, size = 1) +
  labs(
    x = "Mean annual LCS intensity\n(baseline 1981–2010)",
    y = "Mean annual LCS intensity\n(baseline 1991–2020)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 7. Figure 8: histogram panels (c, f, i)
# ------------------------------------------------

# 7.1. Frequency differences (c)
baseline_hist <- baseline_comp %>%
  dplyr::mutate(region = ifelse(lat >= 0,
                                "Northern hemisphere",
                                "Southern hemisphere"))

breaks_freq <- pretty(range(baseline_hist$freq_diff, na.rm = TRUE), n = 60)

hist_n_freq <- hist(
  baseline_hist$freq_diff[baseline_hist$region == "Northern hemisphere"],
  breaks = breaks_freq,
  plot   = FALSE
)
hist_s_freq <- hist(
  baseline_hist$freq_diff[baseline_hist$region == "Southern hemisphere"],
  breaks = breaks_freq,
  plot   = FALSE
)

north_max_freq <- max(hist_n_freq$counts)
south_max_freq <- max(hist_s_freq$counts)
k_freq <- north_max_freq / south_max_freq

df_hist_freq <- bind_rows(
  tibble(
    x      = hist_n_freq$mids,
    count  = hist_n_freq$counts,
    region = "Northern hemisphere"
  ),
  tibble(
    x      = hist_s_freq$mids,
    count  = hist_s_freq$counts * k_freq,
    region = "Southern hemisphere"
  )
)

fig8c <- ggplot(df_hist_freq, aes(x = x, y = count, fill = region)) +
  geom_col(position = "identity", alpha = 0.3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  scale_colour_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_fill_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_y_continuous(
    name = "Number of lakes\n(Northern hemisphere)",
    sec.axis = sec_axis(
      ~ . / k_freq,
      name = "Number of lakes\n(Southern hemisphere)"
    )
  ) +
  labs(
    x = expression(atop(
      Delta * " mean annual LCS frequency",
      "(1991–2020 minus 1981–2010 baseline)"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_text(colour = "purple4"),
    axis.text.y        = element_text(colour = "purple4"),
    axis.title.y.right = element_text(colour = "darkorange2"),
    axis.text.y.right  = element_text(colour = "darkorange2"),
    legend.position    = "none"
  )

# 7.2. Duration differences (f)
baseline_hist <- baseline_comp %>%
  dplyr::mutate(region = ifelse(lat >= 0,
                                "Northern hemisphere",
                                "Southern hemisphere"))

breaks_dur <- pretty(range(baseline_hist$dur_diff, na.rm = TRUE), n = 60)

hist_n_dur <- hist(
  baseline_hist$dur_diff[baseline_hist$region == "Northern hemisphere"],
  breaks = breaks_dur,
  plot   = FALSE
)
hist_s_dur <- hist(
  baseline_hist$dur_diff[baseline_hist$region == "Southern hemisphere"],
  breaks = breaks_dur,
  plot   = FALSE
)

north_max_dur <- max(hist_n_dur$counts)
south_max_dur <- max(hist_s_dur$counts)
k_dur <- north_max_dur / south_max_dur

df_hist_dur <- bind_rows(
  tibble(
    x      = hist_n_dur$mids,
    count  = hist_n_dur$counts,
    region = "Northern hemisphere"
  ),
  tibble(
    x      = hist_s_dur$mids,
    count  = hist_s_dur$counts * k_dur,
    region = "Southern hemisphere"
  )
)

fig8f <- ggplot(df_hist_dur, aes(x = x, y = count, fill = region)) +
  geom_col(position = "identity", alpha = 0.3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  scale_colour_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_fill_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_y_continuous(
    name = "Number of lakes\n(Northern hemisphere)",
    sec.axis = sec_axis(
      ~ . / k_dur,
      name = "Number of lakes\n(Southern hemisphere)"
    )
  ) +
  labs(
    x = expression(atop(
      Delta * " mean annual LCS duration",
      "(1991–2020 minus 1981–2010 baseline)"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_text(colour = "purple4"),
    axis.text.y        = element_text(colour = "purple4"),
    axis.title.y.right = element_text(colour = "darkorange2"),
    axis.text.y.right  = element_text(colour = "darkorange2"),
    legend.position    = "none"
  )

# 7.3. Intensity differences (i)
baseline_hist <- baseline_comp %>%
  dplyr::mutate(region = ifelse(lat >= 0,
                                "Northern hemisphere",
                                "Southern hemisphere"))

breaks_int <- pretty(range(baseline_hist$int_cum_diff, na.rm = TRUE), n = 60)

hist_n_int <- hist(
  baseline_hist$int_cum_diff[baseline_hist$region == "Northern hemisphere"],
  breaks = breaks_int,
  plot   = FALSE
)
hist_s_int <- hist(
  baseline_hist$int_cum_diff[baseline_hist$region == "Southern hemisphere"],
  breaks = breaks_int,
  plot   = FALSE
)

north_max_int <- max(hist_n_int$counts)
south_max_int <- max(hist_s_int$counts)
k_int <- north_max_int / south_max_int

df_hist_int <- bind_rows(
  tibble(
    x      = hist_n_int$mids,
    count  = hist_n_int$counts,
    region = "Northern hemisphere"
  ),
  tibble(
    x      = hist_s_int$mids,
    count  = hist_s_int$counts * k_int,
    region = "Southern hemisphere"
  )
)

fig8i <- ggplot(df_hist_int, aes(x = x, y = count, fill = region)) +
  geom_col(position = "identity", alpha = 0.3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  scale_colour_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_fill_manual(
    values = c("Northern hemisphere" = "purple4",
               "Southern hemisphere" = "darkorange2")
  ) +
  scale_y_continuous(
    name = "Number of lakes\n(Northern hemisphere)",
    sec.axis = sec_axis(
      ~ . / k_int,
      name = "Number of lakes\n(Southern hemisphere)"
    )
  ) +
  labs(
    x = expression(atop(
      Delta * " mean annual LCS intensity",
      "(1991–2020 minus 1981–2010 baseline)"
    ))
  ) +
  theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_text(colour = "purple4"),
    axis.text.y        = element_text(colour = "purple4"),
    axis.title.y.right = element_text(colour = "darkorange2"),
    axis.text.y.right  = element_text(colour = "darkorange2"),
    legend.position    = "none"
  )

# ------------------------------------------------
# 8. Figure 8: global maps (a, d, g)
# ------------------------------------------------

# Palette for clamped differences
cols11      <- hcl.colors(10, "PiYG", rev = TRUE)
cols_custom <- cols11[-c(5, 6)]

# (a) Map of frequency differences
baseline_map_freq <- baseline_comp %>%
  dplyr::mutate(freq_diff_clamped = pmax(pmin(freq_diff, 1), -1))

fig8a <- ggplot(baseline_map_freq, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  geom_point(aes(colour = freq_diff_clamped), size = 0.6) +
  coord_quickmap() +
  scale_colour_stepsn(
    colours = cols_custom,
    limits  = c(-1, 1),
    breaks  = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    name    = expression(
      atop(
        Delta ~ "LCS freq",
        atop("1991–2020 minus", "1981–2010")
      )
    )
  ) +
  scale_y_continuous(
    limits = c(min(baseline_map_freq$lat), max(baseline_map_freq$lat)),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(min(baseline_map_freq$lon), max(baseline_map_freq$lon))
  ) +
  annotate(
    geom  = "text",
    x     = -175,
    y     = 10,
    label = "Equator",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)"
  ) +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank()
  )

# (d) Map of duration differences
baseline_map_dur <- baseline_comp %>%
  dplyr::mutate(dur_diff_clamped = pmax(pmin(dur_diff, 1), -1))

fig8d <- ggplot(baseline_map_dur, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  geom_point(aes(colour = dur_diff_clamped), size = 0.5) +
  coord_quickmap() +
  scale_colour_stepsn(
    colours = cols_custom,
    limits  = c(-1, 1),
    breaks  = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    name    = expression(
      atop(
        Delta ~ "LCS dur",
        atop("1991–2020 minus", "1981–2010")
      )
    )
  ) +
  scale_y_continuous(
    limits = c(min(baseline_map_dur$lat), max(baseline_map_dur$lat)),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(min(baseline_map_dur$lon), max(baseline_map_dur$lon))
  ) +
  annotate(
    geom  = "text",
    x     = -175,
    y     = 10,
    label = "Equator",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)"
  ) +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank()
  )

# (g) Map of intensity differences
baseline_map_int <- baseline_comp %>%
  dplyr::mutate(int_diff_clamped = pmax(pmin(int_cum_diff, 1), -1))

fig8g <- ggplot(baseline_map_int, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  geom_point(aes(colour = int_diff_clamped), size = 0.5) +
  coord_quickmap() +
  scale_colour_stepsn(
    colours = cols_custom,
    limits  = c(-1, 1),
    breaks  = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    name    = expression(
      atop(
        Delta ~ "LCS int",
        atop("1991–2020 minus", "1981–2010")
      )
    )
  ) +
  scale_y_continuous(
    limits = c(min(baseline_map_int$lat), max(baseline_map_int$lat)),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(min(baseline_map_int$lon), max(baseline_map_int$lon))
  ) +
  annotate(
    geom  = "text",
    x     = -175,
    y     = 10,
    label = "Equator",
    colour = "black"
  ) +
  theme_bw() +
  labs(
    x = "Longitude (°)",
    y = "Latitude (°)"
  ) +
  theme(
    panel.background   = element_rect(fill = "aliceblue"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank()
  )

# ------------------------------------------------
# 9. Assemble Figure 8 and save
# ------------------------------------------------

# Row 1: (a) map + (b, c) scatter/hist
row1 <- plot_grid(
  fig8a, fig8b, fig8c,
  ncol       = 3,
  rel_widths = c(2, 1, 1),
  labels     = c("(a)", "(b)", "(c)")
)

# Row 2: (d, e, f)
row2 <- plot_grid(
  fig8d, fig8e, fig8f,
  ncol       = 3,
  rel_widths = c(2, 1, 1),
  labels     = c("(d)", "(e)", "(f)")
)

# Row 3: (g, h, i)
row3 <- plot_grid(
  fig8g, fig8h, fig8i,
  ncol       = 3,
  rel_widths = c(2, 1, 1),
  labels     = c("(g)", "(h)", "(i)")
)

fig8 <- plot_grid(
  row1, row2, row3,
  ncol  = 1,
  align = "v"
)

ggsave(
  filename = file.path(fig_dir, "fig8_baseline_sensitivity.pdf"),
  plot     = fig8,
  device   = pdf,
  width    = 18,
  height   = 9,
  units    = "in",
  dpi      = 300
)