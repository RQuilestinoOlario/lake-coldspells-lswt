#!/usr/bin/env Rscript
# ================================================================
# Example lake cold-spell event (Figure 2)
# ------------------------------------------------
# Description:
#   - Detect lake cold-spells (LCS) for lake 1 using a 1981–2010
#     baseline.
#   - Select a specific event (by event_no) as an illustrative example.
#   - Plot temperature, climatology, threshold, and shaded LCS area.
#
# Requirements:
#   - LSWT data in HDF5/MAT format: "daily_LSWT_data.mat"
#   - Lake and date metadata objects available in the workspace:
#       dateinfo: data.frame with column `date`
#   - This script assumes `dateinfo` was produced by the data-loading
#     script (e.g. 01_load_lswt_metadata.R).
#
# Author: Raven Quilestino-Olario
# Date:   2025-12-10
# License: See repository-level LICENSE file
# ================================================================

suppressPackageStartupMessages({
  library(rhdf5)
  library(heatwaveR)
  library(dplyr)
  library(tibble)
  library(ggplot2)
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

# Full time index (assumes one row per day in dateinfo)
idx_all <- seq_len(nrow(dateinfo))

# ------------------------------------------------
# 2. LCS detection for lake 1
# ------------------------------------------------

# Lake 1, full record in °C
lswt_1 <- h5read(
  file  = f_temp,
  name  = "daily_LSWT_data",
  index = list(1, idx_all)
) - 273.15

ts_lake1 <- tibble::tibble(
  t    = dateinfo$date[idx_all],
  temp = as.numeric(lswt_1)
)

clim_lake1 <- ts2clm(
  data              = ts_lake1,
  x                 = t,
  y                 = temp,
  climatologyPeriod = c("1981-01-01", "2010-12-31"),
  pctile            = 10,
  windowHalfWidth   = 5
)

lcs_lake1 <- detect_event(
  data         = clim_lake1,
  coldSpells   = TRUE,
  min_duration = 5
)

# ------------------------------------------------
# 3. Select an example event and subset climatology
# ------------------------------------------------

# Manually-chosen event number (example)
ev_id <- 49

stopifnot("Requested event_no not found." =
            ev_id %in% lcs_lake1$event$event_no)

ev <- lcs_lake1$event %>%
  dplyr::filter(event_no == ev_id)

# Subset the climatology around that event
clim_sub <- lcs_lake1$clim %>%
  dplyr::filter(
    t >= ev$date_start - 10,
    t <= ev$date_end   + 10
  ) %>%
  dplyr::mutate(
    in_event = t >= ev$date_start & t <= ev$date_end
  )

max_y <- max(clim_sub$temp, clim_sub$thresh, na.rm = TRUE)
min_y <- min(clim_sub$temp, clim_sub$thresh, na.rm = TRUE)

# ------------------------------------------------
# 4. Plot: temperature, climatology, threshold, shaded LCS
# ------------------------------------------------

fig2 <- ggplot(clim_sub, aes(x = t)) +
  # Shaded cold-spell area (only during the event)
  geom_ribbon(
    data  = dplyr::filter(clim_sub, in_event),
    aes(ymin = temp, ymax = thresh),
    fill  = "lightskyblue2",
    alpha = 0.8
  ) +
  # Lines
  geom_line(aes(y = seas,    colour = "Climatology"), linewidth = 1.0) +
  geom_line(aes(y = thresh,  colour = "Threshold"),   linewidth = 1.0) +
  geom_line(aes(y = temp,    colour = "Temperature"), linewidth = 0.7) +
  scale_colour_manual(
    name   = NULL,
    values = c(
      "Temperature" = "black",
      "Climatology" = "darkorange1",
      "Threshold"   = "deepskyblue4"
    )
  ) +
  scale_x_date(
    expand      = c(0, 0),
    date_labels = "%b %Y"
  ) +
  scale_y_continuous(
    limits = c(min_y - 0.5, max_y + 2),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Lake surface water temperature (°C)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  )

# ------------------------------------------------
# 5. Save Figure 2
# ------------------------------------------------

ggsave(
  filename = file.path(fig_dir, "fig2_example_lcs_event.png"),
  plot     = fig2,
  device   = png,
  width    = 6000,
  height   = 4012.5,
  units    = "px",
  dpi      = 1200
)