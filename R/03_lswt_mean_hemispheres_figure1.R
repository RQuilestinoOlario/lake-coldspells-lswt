#!/usr/bin/env Rscript
# ================================================================
# Mean LSWT and hemispheric temperature structure (Figure 1)
# ------------------------------------------------
# Description:
#   - Compute mean lake surface water temperature (LSWT) over 1981–2020
#     for all lakes.
#   - Map global mean LSWT (panel a).
#   - Show meridional structure and hemispheric division (panel b).
#
# Requirements:
#   - LSWT data in HDF5/MAT format: "daily_LSWT_data.mat"
#   - Lake and date metadata objects available in the workspace:
#       lakeinfo: data.frame with at least columns Hylak_id, lat, lon
#       dateinfo: data.frame with one row per daily time step
#     (e.g. loaded from RDS files before running this script).
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(maps)
})

# ------------------------------------------------
# 1. Paths and basic setup
# ------------------------------------------------

# Directory containing the LSWT input files
data_dir <- "your/path/to/data"

# LSWT data file (daily LSWT, K)
f_temp <- file.path(data_dir, "daily_LSWT_data.mat")

stopifnot("Temperature file not found." = file.exists(f_temp))

# Figure output directory
fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. Indices and chunking
# ------------------------------------------------

# Indices for climatology period
# Here we use the full 1981–2020 record (one row per day in dateinfo)
idx_clim <- seq_len(nrow(dateinfo))

n_lakes    <- nrow(lakeinfo)
chunk_size <- 1000
lake_ids   <- seq_len(n_lakes)
chunks     <- split(lake_ids, ceiling(seq_along(lake_ids) / chunk_size))

# ------------------------------------------------
# 3. Compute mean LSWT (1981–2020) per lake
# ------------------------------------------------

lswt_mean_full <- numeric(n_lakes)

for (ch in chunks) {
  message("Chunk ", min(ch), "–", max(ch))
  
  tmp <- h5read(
    file  = f_temp,
    name  = "daily_LSWT_data",
    index = list(ch, idx_clim)
  ) - 273.15  # K -> °C
  
  lswt_mean_full[ch] <- rowMeans(tmp, na.rm = TRUE)
}

lakeinfo_temp <- lakeinfo %>%
  dplyr::mutate(lswt_mean_8120 = lswt_mean_full)

# ------------------------------------------------
# 4. Base map and palette
# ------------------------------------------------

map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>%
  dplyr::rename(lon = long)

cols10       <- hcl.colors(10, "Cividis")
cols_custom2 <- cols10

# ------------------------------------------------
# 5. Panel (a): global map of mean LSWT
# ------------------------------------------------

fig1a <- ggplot(lakeinfo_temp, aes(x = lon, y = lat)) +
  geom_polygon(
    data   = map_base,
    aes(group = group),
    colour = NA,
    fill   = "grey90"
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  annotate(
    geom  = "text",
    x     = -175,
    y     = 10,
    label = "Equator",
    colour = "black"
  ) +
  geom_point(
    aes(colour = lswt_mean_8120),
    size  = 0.35,
    alpha = 0.8
  ) +
  coord_quickmap() +
  scale_x_continuous(limits = range(lakeinfo_temp$lon, na.rm = TRUE)) +
  scale_y_continuous(
    limits = range(lakeinfo_temp$lat, na.rm = TRUE),
    expand = c(0, 0)
  ) +
  scale_colour_gradientn(
    colours = cols_custom2,
    name    = "Mean LSWT\n1981–2020\n(°C)"
  ) +
  theme_bw() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------
# 6. Panel (b): latitude vs mean LSWT by hemisphere
# ------------------------------------------------

lakeinfo_temp2 <- lakeinfo_temp %>%
  dplyr::mutate(
    hemisphere = ifelse(lat >= 0,
                        "Northern hemisphere",
                        "Southern hemisphere")
  )

fig1b <- ggplot(lakeinfo_temp2,
                aes(x = lat, y = lswt_mean_8120)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "black") +
  annotate(
    geom  = "text",
    x     = -2.5,
    y     = 5,
    label = "Equator",
    colour = "black",
    angle = 90
  ) +
  geom_point(
    aes(colour = hemisphere),
    alpha = 0.5,
    size  = 0.4
  ) +
  geom_smooth(se = FALSE, colour = "black") +
  scale_colour_manual(
    values = c(
      "Northern hemisphere" = "purple4",
      "Southern hemisphere" = "darkorange2"
    ),
    name = NULL
  ) +
  theme_bw() +
  labs(
    x = "Latitude (°)",
    y = "Mean LSWT (°C, 1981–2020)"
  ) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = c(0.85, 0.8),
    legend.background = element_rect(
      fill   = "white",
      colour = "black",
      size   = 0.3
    )
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

# ------------------------------------------------
# 7. Combine panels and save Figure 1
# ------------------------------------------------

fig1 <- plot_grid(
  fig1a, fig1b,
  labels = c("(a)", "(b)"),
  ncol   = 1,
  nrow   = 2,
  align  = "v"
)

theme_set(theme_minimal(base_size = 10))

ggsave(
  filename = file.path(fig_dir, "fig1_lswt_map_latitude.pdf"),
  plot     = fig1,
  device   = pdf,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)