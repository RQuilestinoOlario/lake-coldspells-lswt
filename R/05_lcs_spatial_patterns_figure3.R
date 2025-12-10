#!/usr/bin/env Rscript
# ================================================================
# Spatial patterns of LCS metrics (Figure 3)
# ------------------------------------------------
# Description:
#   - Using lake cold-spell (LCS) metrics from the 1981–2010 baseline
#     (lcs_8110_geo), compute mean values over 1981–2020:
#       * Mean annual LCS frequency (events/year)
#       * Mean LCS duration (days/event)
#       * Mean LCS intensity (°C/event)
#   - Map spatial patterns of these metrics (panels a–c).
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
# 1. Figure/auxiliary paths (optional)
# ------------------------------------------------

fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. Mean LCS metrics over 1981–2020
# ------------------------------------------------

lcs_clim_8120 <- lcs_8110_geo %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::group_by(Hylak_id, lake_idx, lat, lon) %>%
  dplyr::summarise(
    freq_mean = mean(n_events,            na.rm = TRUE),  # events/year
    dur_mean  = mean(duration_mean,       na.rm = TRUE),  # days/event
    int_mean  = mean(intensity_mean_mean, na.rm = TRUE),  # °C/event
    .groups   = "drop"
  )

# ------------------------------------------------
# 3. Base map and colour palette
# ------------------------------------------------

map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long)

cols10       <- hcl.colors(10, "Cividis")
cols_custom2 <- cols10

# ------------------------------------------------
# 4. Panel (a): mean annual LCS frequency
# ------------------------------------------------

fig3a <- ggplot(lcs_clim_8120, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed", colour = "black") +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = freq_mean),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_clim_8120$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_clim_8120$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = cols_custom2,
    name    = "Mean annual\nLCS frequency\n1981–2020\n(events/year)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 5. Panel (b): mean LCS duration (days/event)
# ------------------------------------------------

fig3b <- ggplot(lcs_clim_8120, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed", colour = "black") +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = dur_mean),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_clim_8120$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_clim_8120$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = cols_custom2,
    name    = "Mean LCS\nduration\n1981–2020\n(days/event)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 6. Panel (c): mean LCS intensity (°C/event)
# ------------------------------------------------

fig3c <- ggplot(lcs_clim_8120, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0),
             linetype = "dashed", colour = "black") +
  annotate("text", x = -175, y = 10, label = "Equator", colour = "black") +
  geom_point(
    aes(colour = int_mean),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lcs_clim_8120$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lcs_clim_8120$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = rev(cols_custom2),
    name    = "Mean LCS\nintensity\n1981–2020\n(°C/event)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 7. Combine panels and save Figure 3
# ------------------------------------------------

fig3 <- plot_grid(
  fig3a, fig3b, fig3c,
  labels = c("(a)", "(b)", "(c)"),
  ncol   = 1,
  nrow   = 3,
  align  = "v"
)

theme_set(theme_minimal(base_size = 10))

ggsave(
  filename = file.path(fig_dir, "fig3_annual_spatial_patterns.pdf"),
  plot     = fig3,
  device   = pdf,
  width    = 8,
  height   = 9,
  units    = "in",
  dpi      = 300
)