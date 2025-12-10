#!/usr/bin/env Rscript
# ================================================================
# Latitudinal and temperature-band structure of LCS metrics (Figure 5)
# ------------------------------------------------
# Description:
#   - Bin lakes by latitude and compute band-mean LCS metrics:
#       * Mean annual LCS frequency (events/year)
#       * Mean LCS duration (days/event)
#       * Mean LCS intensity (°C/event)
#   - Bin lakes by mean LSWT (1981–2020) and summarise the same
#     metrics with boxplots across temperature bands.
#
# Requirements:
#   - `lcs_clim_8120` from 05_lcs_spatial_patterns_figure3.R, with:
#       Hylak_id, lat, lon, freq_mean, dur_mean, int_mean
#   - `lakeinfo_temp` from 03_lswt_mean_hemispheres_figure1.R, with:
#       Hylak_id, lswt_mean_8120
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# ------------------------------------------------
# 1. Figure directory
# ------------------------------------------------

fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. Latitude bands: zonal means of LCS metrics
# ------------------------------------------------

lcs_lat <- lcs_clim_8120 %>%
  dplyr::mutate(
    lat_band = cut(
      lat,
      breaks = seq(-60, 80, by = 5),  # adjust range/bandwidth as needed
      include.lowest = TRUE
    )
  ) %>%
  dplyr::group_by(lat_band) %>%
  dplyr::summarise(
    lat_mid   = mean(lat,       na.rm = TRUE),
    freq_mean = mean(freq_mean, na.rm = TRUE),
    dur_mean  = mean(dur_mean,  na.rm = TRUE),
    int_mean  = mean(int_mean,  na.rm = TRUE),
    n_lakes   = dplyr::n(),
    .groups   = "drop"
  )

# Panel (a): frequency vs latitude
fig5a <- ggplot(lcs_lat, aes(x = lat_mid, y = freq_mean)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red") +
  geom_line() +
  geom_point() +
  labs(
    x = "Latitude (°)",
    y = "Mean annual LCS\nfrequency (events/year)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel (b): duration vs latitude
fig5b <- ggplot(lcs_lat, aes(x = lat_mid, y = dur_mean)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red") +
  annotate(
    geom  = "text",
    x     = -2.5,
    y     = 15,
    label = "Equator",
    colour = "red",
    angle = 90
  ) +
  geom_line() +
  geom_point() +
  labs(
    x = "Latitude (°)",
    y = "Mean LCS duration\n(days/event)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel (c): intensity vs latitude
fig5c <- ggplot(lcs_lat, aes(x = lat_mid, y = int_mean)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red") +
  geom_line() +
  geom_point() +
  labs(
    x = "Latitude (°)",
    y = "Mean LCS intensity\n(°C/event)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 3. Temperature bands: boxplots vs mean LSWT
# ------------------------------------------------

lcs_lswt_band <- lcs_clim_8120 %>%
  dplyr::left_join(
    lakeinfo_temp %>% dplyr::select(Hylak_id, lswt_mean_8120),
    by = "Hylak_id"
  ) %>%
  dplyr::mutate(
    temp_band = cut(
      lswt_mean_8120,
      breaks = c(-Inf, 5, 10, 15, 20, 25, Inf),
      labels = c("<5 °C", "5–10 °C", "10–15 °C",
                 "15–20 °C", "20–25 °C", ">25 °C"),
      right  = TRUE
    )
  )

# Panel (d): frequency vs temp band
fig5d <- ggplot(lcs_lswt_band, aes(x = temp_band, y = freq_mean)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(
    x = "Mean LSWT band (1981–2020)",
    y = "Mean annual LCS\nfrequency (events/year)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel (e): duration vs temp band
fig5e <- ggplot(lcs_lswt_band, aes(x = temp_band, y = dur_mean)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(
    x = "Mean LSWT band (1981–2020)",
    y = "Mean LCS duration\n(days/event)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel (f): intensity vs temp band
fig5f <- ggplot(lcs_lswt_band, aes(x = temp_band, y = int_mean)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(
    x = "Mean LSWT band (1981–2020)",
    y = "Mean LCS intensity\n(°C/event)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 4. Combine panels and save Figure 5
# ------------------------------------------------

fig5 <- plot_grid(
  fig5a, fig5d,
  fig5b, fig5e,
  fig5c, fig5f,
  ncol   = 2,
  align  = "hv",
  labels = c("(a)", "(d)", "(b)", "(e)", "(c)", "(f)")
)

save_plot(
  filename    = file.path(fig_dir, "fig5_climateboxplots.pdf"),
  plot        = fig5,
  base_width  = 12,
  base_height = 8
)