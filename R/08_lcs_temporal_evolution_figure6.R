#!/usr/bin/env Rscript
# ================================================================
# Temporal evolution of LCS metrics by latitude band (Figure 6)
# ------------------------------------------------
# Description:
#   - Using lake cold-spell (LCS) metrics from the 1981–2010 baseline
#     (lcs_8110_geo), compute annual mean LCS metrics by latitude band:
#       * Mean annual LCS frequency (events/year)
#       * Mean LCS duration (days/event)
#       * Mean LCS intensity (°C/event)
#   - Also compute global annual means.
#   - Plot temporal evolution per band with:
#       - line: frequency
#       - points: duration (scaled to same axis) with size/colour
#         indicating intensity.
#
# Requirements:
#   - Object `lcs_8110_geo` in memory, with at least:
#       year, lat,
#       n_events, duration_mean, intensity_mean_mean.
#   - Typically produced in 02_lcs_baseline_sensitivity_figure8.R.
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# ------------------------------------------------
# 1. Figure directory
# ------------------------------------------------

fig_dir <- "your/path/to/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# 2. Latitude-band annual means
# ------------------------------------------------

lcs_band_year <- lcs_8110_geo %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::mutate(
    lat_band = dplyr::case_when(
      lat >= 60        ~ "60–90°N",
      lat >= 30        ~ "30–60°N",
      lat >= 0         ~ "0–30°N",
      lat >  -30       ~ "0–30°S",
      lat >  -60       ~ "30–60°S",
      TRUE             ~ "60–90°S"
    )
  ) %>%
  dplyr::group_by(lat_band, year) %>%
  dplyr::summarise(
    freq_mean = mean(n_events,            na.rm = TRUE),
    dur_mean  = mean(duration_mean,       na.rm = TRUE),
    int_mean  = mean(intensity_mean_mean, na.rm = TRUE),
    .groups   = "drop"
  )

# ------------------------------------------------
# 3. Global annual means
# ------------------------------------------------

lcs_global_year <- lcs_8110_geo %>%
  dplyr::filter(year >= 1981, year <= 2020) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    freq_mean = mean(n_events,            na.rm = TRUE),
    dur_mean  = mean(duration_mean,       na.rm = TRUE),
    int_mean  = mean(intensity_mean_mean, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(lat_band = "Global")

# ------------------------------------------------
# 4. Combine bands and set facet order (drop empty 60–90°S)
# ------------------------------------------------

lcs_band_year_all <- dplyr::bind_rows(lcs_band_year, lcs_global_year) %>%
  dplyr::filter(lat_band != "60–90°S") %>%  # drop empty Antarctic band
  dplyr::mutate(
    lat_band = factor(
      lat_band,
      levels = c("60–90°N", "30–60°N", "0–30°N",
                 "0–30°S", "30–60°S", "Global")
    )
  )

# ------------------------------------------------
# 5. Global ranges and transforms for dual y-axis
# ------------------------------------------------

freq_rng <- range(lcs_band_year_all$freq_mean, na.rm = TRUE)
dur_rng  <- range(lcs_band_year_all$dur_mean,  na.rm = TRUE)

dur_to_freq <- function(d) {
  (d - dur_rng[1]) / diff(dur_rng) * diff(freq_rng) + freq_rng[1]
}

freq_to_dur <- function(y) {
  (y - freq_rng[1]) / diff(freq_rng) * diff(dur_rng) + dur_rng[1]
}

lcs_band_year_all <- lcs_band_year_all %>%
  dplyr::mutate(dur_scaled = dur_to_freq(dur_mean))

# ------------------------------------------------
# 6. Intensity bins and colours (Okabe–Ito)
# ------------------------------------------------

okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7"  # reddish purple
)

lcs_band_year_all <- lcs_band_year_all %>%
  dplyr::mutate(
    int_abs = abs(int_mean),
    int_bin = cut(
      int_abs,
      breaks = c(0, 2, 2.5, 3, 3.5, 4, Inf),
      labels = c("≤2", "2–2.5", "2.5–3", "3–3.5", "3.5–4", ">4"),
      right  = TRUE
    )
  )

# ------------------------------------------------
# 7. Figure 6: temporal evolution by latitude band
# ------------------------------------------------

fig6 <- ggplot(lcs_band_year_all, aes(x = year)) +
  # line: mean annual frequency
  geom_line(
    aes(y = freq_mean),
    colour   = "black",
    linewidth = 0.5
  ) +
  # points: duration (scaled) coloured/sized by intensity
  geom_point(
    aes(
      y      = dur_scaled,
      size   = int_bin,
      colour = int_bin
    ),
    alpha = 0.9
  ) +
  facet_wrap(~ lat_band, ncol = 2) +
  scale_y_continuous(
    name    = "Mean annual LCS frequency (events/year)",
    sec.axis = sec_axis(
      ~ freq_to_dur(.),
      name = "Mean LCS duration (days/event)"
    )
  ) +
  scale_size_manual(
    name   = "Mean LCS intensity\n(|°C| per event)",
    values = c(
      "≤2"    = 1.8,
      "2–2.5" = 2.2,
      "2.5–3" = 2.6,
      "3–3.5" = 3.0,
      "3.5–4" = 3.4,
      ">4"    = 3.8
    )
  ) +
  scale_colour_manual(
    name   = "Mean LCS intensity\n(|°C| per event)",
    values = c(
      "≤2"    = okabe_ito[1],
      "2–2.5" = okabe_ito[2],
      "2.5–3" = okabe_ito[3],
      "3–3.5" = okabe_ito[4],
      "3.5–4" = okabe_ito[5],
      ">4"    = okabe_ito[6]
    )
  ) +
  theme_bw() +
  labs(x = "Year") +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "white"),
    legend.position    = "right",
    axis.title.y.right = element_text(colour = "#D55E00"),
    axis.text.y.right  = element_text(colour = "#D55E00"),
    axis.line.y.right  = element_line(colour = "#D55E00"),
    axis.ticks.y.right = element_line(colour = "#D55E00")
  ) +
  guides(
    size   = guide_legend(override.aes = list(shape = 16)),
    colour = "legend"
  )

# ------------------------------------------------
# 8. Save Figure 6
# ------------------------------------------------

ggsave(
  filename = file.path(fig_dir, "fig6_timeevolution.pdf"),
  plot     = fig6,
  device   = pdf,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 300
)