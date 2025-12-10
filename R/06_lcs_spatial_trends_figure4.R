#!/usr/bin/env Rscript
# ================================================================
# Spatial trends in LCS metrics (Figure 4)
# ------------------------------------------------
# Description:
#   - Using lake cold-spell (LCS) metrics from the 1981–2010 baseline
#     (lcs_8110_geo), compute linear trends over 1981–2020 for:
#       * Mean annual LCS frequency (events/decade)
#       * Mean LCS duration (days/event/decade)
#       * Mean LCS intensity (°C/event/decade)
#   - Trends are estimated per lake via simple linear regression
#     against calendar year, scaled to per-decade rates.
#   - Only lakes with at least 20 years of data are retained.
#   - Map spatial patterns of these trends (panels a–c).
#
# Requirements:
#   - Aggregated LCS output with geometry, e.g. from
#       02_lcs_baseline_sensitivity_figure8.R
#     stored as an object `lcs_8110_geo` with at least:
#       Hylak_id, lake_idx, lat, lon, year,
#       n_events, duration_mean, intensity_mean_mean.
#
#   Example of loading from disk:
#       lcs_8110_geo <- readRDS("your/path/to/lcs_8110_geo.rds")
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(maps)
})

# ------------------------------------------------
# 1. Figure directory
# ------------------------------------------------

fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. Restrict to analysis period and compute per-lake trends
# ------------------------------------------------

lcs_trend_input <- lcs_8110_geo %>%
  dplyr::filter(year >= 1981, year <= 2020)

# Per-lake trends; require at least 20 distinct years of data
lcs_trends <- lcs_trend_input %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::filter(dplyr::n_distinct(year) >= 20) %>%
  dplyr::summarise(
    freq_trend_decade = coef(lm(n_events            ~ year))[["year"]] * 10,
    dur_trend_decade  = coef(lm(duration_mean       ~ year))[["year"]] * 10,
    int_trend_decade  = coef(lm(intensity_mean_mean ~ year))[["year"]] * 10,
    n_years           = dplyr::n(),
    .groups           = "drop"
  )

# ------------------------------------------------
# 3. Trend limits and colour palette
# ------------------------------------------------

lim_freq <- quantile(abs(lcs_trends$freq_trend_decade), 0.99, na.rm = TRUE)
lim_dur  <- quantile(abs(lcs_trends$dur_trend_decade),  0.99, na.rm = TRUE)
lim_int  <- quantile(abs(lcs_trends$int_trend_decade),  0.99, na.rm = TRUE)

cols11      <- hcl.colors(10, "PuOr", rev = TRUE)
cols_custom <- cols11[-c(5, 6)]

# Base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long)

# ------------------------------------------------
# 4. Panel (a): Trend in mean annual LCS frequency
# ------------------------------------------------

fig4a <- ggplot(lcs_trends, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed",
    colour   = "black"
  ) +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = freq_trend_decade),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_trends$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_trends$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_custom),
    limits  = c(-lim_freq, lim_freq),
    name    = "Trend in mean\nannual LCS freq\n1981–2020\n(events/decade)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 5. Panel (b): Trend in mean LCS duration
# ------------------------------------------------

fig4b <- ggplot(lcs_trends, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed",
    colour   = "black"
  ) +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = dur_trend_decade),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_trends$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_trends$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_custom),
    limits  = c(-lim_dur, lim_dur),
    name    = "Trend in mean\nLCS duration\n1981–2020\n(days/event/decade)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 6. Panel (c): Trend in mean LCS intensity
# ------------------------------------------------

fig4c <- ggplot(lcs_trends, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed",
    colour   = "black"
  ) +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = int_trend_decade),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_trends$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_trends$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = cols_custom,
    limits  = c(-lim_int, lim_int),
    name    = "Trend in mean\nLCS intensity\n1981–2020\n(°C/event/decade)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 7. Combine panels and save Figure 4
# ------------------------------------------------

fig4 <- plot_grid(
  fig4a, fig4b, fig4c,
  ncol   = 1,
  labels = c("(a)", "(b)", "(c)"),
  align  = "v"
)

theme_set(theme_minimal(base_size = 10))

ggsave(
  filename = file.path(fig_dir, "fig4_trends_spatial_patterns.pdf"),
  plot     = fig4,
  device   = pdf,
  width    = 8,
  height   = 9,
  units    = "in",
  dpi      = 300
)