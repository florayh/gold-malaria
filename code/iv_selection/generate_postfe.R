#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Generate post-FE residualized panels for Stata rlasso IV selection.
#          Shows reviewers how the pre-computed postFE files in data/postfe/
#          were created. OPTIONAL — pre-computed versions are included.
# Author: Flora He (assisted by Claude Code)
# Start date: 2026-06-11
# Code review: NOT REVIEWED
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Configuration -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Input: analysis panels (already in reproduction package)
panels_dir <- "data/panels/"
spill_dir  <- "data/spillover/"

# Output: post-FE filtered panels
out_dir <- "data/postfe/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gp_suffix <- "_2yrgp"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: quantile-based post-FE filtering -----
# 1. Extract interaction columns (suffix _2yrgp)
# 2. Regress on unit FE: lm(x_mat ~ as.factor(unit_id))$residuals
# 3. Compute variance: apply(xr, 2, var)
# 4. Filter: keep columns with variance > quantile threshold
# 5. Return data frame: unit_id + year + treatment + passing IVs + extras
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postfe_quantile <- function(df, unit_var, treatment_var, quantile_val,
                            extra_vars = character(0),
                            iv_pattern = "_2yrgp$") {
  # Identify candidate IVs
  iv_cols <- grep(iv_pattern, names(df), value = TRUE)
  cat(sprintf("  Candidates: %d interaction columns\n", length(iv_cols)))

  # Partial out unit FE
  x_mat <- as.matrix(df[, iv_cols])
  f <- df[[unit_var]]
  xr <- lm(x_mat ~ as.factor(f))$residuals
  xr_var <- apply(xr, 2, var)

  # Apply quantile threshold
  threshold <- quantile(xr_var, quantile_val)
  selected <- names(xr_var[xr_var > threshold])
  cat(sprintf("  Threshold (q=%.2f): %.1f | Passing: %d of %d\n",
              quantile_val, threshold, length(selected), length(xr_var)))

  # Build output: identifiers + treatment + selected IVs + extras
  keep_cols <- unique(c(unit_var, "year", treatment_var, extra_vars, selected))
  keep_cols <- intersect(keep_cols, names(df))
  df[, keep_cols]
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: drop-zeros + cap post-FE filtering -----
# For band/spillover panels where many IVs have zero post-FE variance.
# 1. Drop IVs with zero post-FE variance
# 2. Sort remaining by variance (descending)
# 3. Keep top max_iv IVs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postfe_drop_zeros_cap <- function(df, unit_var, treatment_var, max_iv,
                                  extra_vars = character(0),
                                  iv_pattern = "_2yrgp$") {
  iv_cols <- grep(iv_pattern, names(df), value = TRUE)
  cat(sprintf("  Candidates: %d interaction columns\n", length(iv_cols)))

  x_mat <- as.matrix(df[, iv_cols])
  f <- df[[unit_var]]
  xr <- lm(x_mat ~ as.factor(f))$residuals
  xr_var <- apply(xr, 2, var)

  # Drop zeros, sort by variance, cap
  nonzero <- xr_var[xr_var > 0]
  nonzero <- sort(nonzero, decreasing = TRUE)
  selected <- names(nonzero)[seq_len(min(max_iv, length(nonzero)))]
  cat(sprintf("  Non-zero: %d | After cap %d: %d\n",
              length(nonzero), max_iv, length(selected)))

  keep_cols <- unique(c(unit_var, "year", treatment_var, extra_vars, selected))
  keep_cols <- intersect(keep_cols, names(df))
  df[, keep_cols]
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Group A: Main 03-19 panels -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("=== Group A: Main 03-19 panels ===\n")

# A1: Municipality gold+neighbor (q=0.25)
cat("\n--- A1: Municipality gold+neighbor ---\n")
mu_df <- read_csv(paste0(panels_dir, "mu_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
mu_postfe <- postfe_quantile(mu_df, "cd_mun", "gold_mining_area", 0.25,
                             extra_vars = c("pop", "malaria_allpop_api",
                                            "forest_area_km", "agriculture_area_km",
                                            "gdp_per_capita", "mean_temp", "annual_precip"))
write_csv(mu_postfe, paste0(out_dir, "mu_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_postFE.csv"))
cat(sprintf("  Saved: %d rows x %d cols\n", nrow(mu_postfe), ncol(mu_postfe)))

# A2: Municipality drop Itaituba + Jacareacanga (q=0.25)
cat("\n--- A2: Municipality drop Itaituba+Jacareacanga ---\n")
mu_dropij <- mu_df %>% filter(!cd_mun %in% c(150360, 150375))
mu_dropij_postfe <- postfe_quantile(mu_dropij, "cd_mun", "gold_mining_area", 0.25,
                                    extra_vars = c("pop", "malaria_allpop_api",
                                                   "forest_area_km", "agriculture_area_km",
                                                   "gdp_per_capita", "mean_temp", "annual_precip"))
write_csv(mu_dropij_postfe, paste0(out_dir, "mu_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_dropIJ_postFE.csv"))
cat(sprintf("  Saved: %d rows x %d cols\n", nrow(mu_dropij_postfe), ncol(mu_dropij_postfe)))

# A3: DSEI gold+neighbor (q=0.45)
cat("\n--- A3: DSEI gold+neighbor ---\n")
dsei_df <- read_csv(paste0(panels_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
dsei_postfe <- postfe_quantile(dsei_df, "seqid", "goldmine_area", 0.45,
                               extra_vars = c("population_polobase", "malaria_api_by_polobase",
                                              "mean_temp", "annual_total_precip", "forest"))
write_csv(dsei_postfe, paste0(out_dir, "dsei_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_postFE.csv"))
cat(sprintf("  Saved: %d rows x %d cols\n", nrow(dsei_postfe), ncol(dsei_postfe)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Group B: Extended panels (03-24) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Group B: Extended panels (03-24) ===\n")

# B1: Municipality extended (q=0.25)
cat("\n--- B1: Municipality extended ---\n")
mu_ext <- read_csv(paste0(panels_dir, "mu_panel_goldneighbor_full.csv"), show_col_types = FALSE)
mu_ext_postfe <- postfe_quantile(mu_ext, "cd_mun", "gold_mining_area", 0.25,
                                 extra_vars = c("pop", "malaria_allpop_api",
                                                "forest_area_km", "agriculture_area_km",
                                                "gdp_per_capita", "mean_temp", "annual_precip"))
write_csv(mu_ext_postfe, paste0(out_dir, "mu_panel_2yrgp_rockarea_20260316_goldneighbor_postFE.csv"))
cat(sprintf("  Saved: %d rows x %d cols\n", nrow(mu_ext_postfe), ncol(mu_ext_postfe)))

# B2: DSEI extended (q=0.45)
cat("\n--- B2: DSEI extended ---\n")
dsei_ext <- read_csv(paste0(panels_dir, "dsei_panel_goldneighbor_full.csv"), show_col_types = FALSE)
dsei_ext_postfe <- postfe_quantile(dsei_ext, "seqid", "goldmine_area", 0.45,
                                   extra_vars = c("population_polobase", "malaria_api_by_polobase",
                                                  "mean_temp", "annual_total_precip", "forest"))
write_csv(dsei_ext_postfe, paste0(out_dir, "dsei_panel_2yrgp_rockarea_20260316_goldneighbor_postFE.csv"))
cat(sprintf("  Saved: %d rows x %d cols\n", nrow(dsei_ext_postfe), ncol(dsei_ext_postfe)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Group C: Distance bands x7 (drop zeros + cap 80) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Group C: Distance band panels ===\n")

band_names <- c("0_50km", "50_100km", "100_150km", "150_200km",
                "200_250km", "250_300km", "300_350km")

for (band in band_names) {
  cat(sprintf("\n--- Band: %s ---\n", band))
  src_file <- paste0(spill_dir, "mu_panel_band_", band, ".csv")
  if (!file.exists(src_file)) {
    cat("  SKIP (not found):", src_file, "\n")
    next
  }
  band_df <- read_csv(src_file, show_col_types = FALSE)
  treatment_var <- paste0("gold_mining_area_nb", sub("km$", "", band), "km")
  band_postfe <- postfe_drop_zeros_cap(band_df, "cd_mun", treatment_var, 80)
  write_csv(band_postfe, paste0(out_dir, "mu_panel_band_", band, "_postFE.csv"))
  cat(sprintf("  Saved: %d rows x %d cols\n", nrow(band_postfe), ncol(band_postfe)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Group D: Clean control + direct spillover (q=0.25) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Group D: Clean control + direct spillover ===\n")

# D1: Clean control (q=0.25)
cat("\n--- D1: Clean control ---\n")
cc_file <- paste0(spill_dir, "mu_panel_clean_control_300km.csv")
if (file.exists(cc_file)) {
  cc_df <- read_csv(cc_file, show_col_types = FALSE)
  cc_postfe <- postfe_quantile(cc_df, "cd_mun", "gold_mining_area", 0.25,
                               extra_vars = c("pop", "malaria_allpop_api",
                                              "forest_area_km", "agriculture_area_km",
                                              "gdp_per_capita", "mean_temp", "annual_precip"))
  write_csv(cc_postfe, paste0(out_dir, "mu_panel_clean_control_300km_postFE.csv"))
  cat(sprintf("  Saved: %d rows x %d cols\n", nrow(cc_postfe), ncol(cc_postfe)))
}

# D2: Direct spillover (q=0.25, applied separately to own and neighbor IV pools)
cat("\n--- D2: Direct spillover ---\n")
ds_file <- paste0(spill_dir, "mu_panel_gold_direct_spillover_200km.csv")
if (file.exists(ds_file)) {
  ds_df <- read_csv(ds_file, show_col_types = FALSE)

  # Separate own vs neighbor IVs
  all_ivs <- grep("_2yrgp$", names(ds_df), value = TRUE)
  nb_ivs  <- grep("_nb0_200km_2yrgp$", all_ivs, value = TRUE)
  own_ivs <- setdiff(all_ivs, nb_ivs)
  cat(sprintf("  Total IVs: %d (own: %d, neighbor: %d)\n",
              length(all_ivs), length(own_ivs), length(nb_ivs)))

  # Filter own pool
  x_own <- as.matrix(ds_df[, own_ivs])
  f <- ds_df$cd_mun
  xr_own <- lm(x_own ~ as.factor(f))$residuals
  var_own <- apply(xr_own, 2, var)
  thresh_own <- quantile(var_own, 0.25)
  sel_own <- names(var_own[var_own > thresh_own])

  # Filter neighbor pool
  x_nb <- as.matrix(ds_df[, nb_ivs])
  xr_nb <- lm(x_nb ~ as.factor(f))$residuals
  var_nb <- apply(xr_nb, 2, var)
  thresh_nb <- quantile(var_nb, 0.25)
  sel_nb <- names(var_nb[var_nb > thresh_nb])

  cat(sprintf("  Own passing: %d of %d | Neighbor passing: %d of %d\n",
              length(sel_own), length(own_ivs), length(sel_nb), length(nb_ivs)))

  extra_vars <- c("pop", "malaria_allpop_api",
                  "forest_area_km", "agriculture_area_km",
                  "gdp_per_capita", "mean_temp", "annual_precip")
  keep_cols <- unique(c("cd_mun", "year", "gold_mining_area",
                        "gold_mining_area_nb0_200km", extra_vars,
                        sel_own, sel_nb))
  keep_cols <- intersect(keep_cols, names(ds_df))
  ds_postfe <- ds_df[, keep_cols]
  write_csv(ds_postfe, paste0(out_dir, "mu_panel_gold_direct_spillover_200km_postFE.csv"))
  cat(sprintf("  Saved: %d rows x %d cols\n", nrow(ds_postfe), ncol(ds_postfe)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Group E: DSEI falsification subsamples (q=0.45) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Group E: DSEI falsification subsamples ===\n")

# Merge DSEI panel with health panel
dsei_main <- read_csv(paste0(panels_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
hp <- read_csv(paste0(panels_dir, "dsei_health_panel.csv"), show_col_types = FALSE)
hp$polo_base_lower <- tolower(hp$polo_base)
hp$dsei_lower <- tolower(hp$dsei)

# Merge on (year, dsei, polo) to avoid duplicates from same polo name in different DSEIs
merged <- merge(dsei_main, hp,
                by.x = c("year", "dsei_gesta", "polo"),
                by.y = c("year", "dsei_lower", "polo_base_lower"),
                all.x = TRUE)
stopifnot(nrow(merged) == nrow(dsei_main))

# Define subsample: births >= 5 with non-missing IMR (matching 03_analysis_dsei.R)
subsamples <- list(
  mortality = merged %>% filter(!is.na(birth_total_live_births) & birth_total_live_births >= 5 &
                                  !is.na(child_mort_imr))
)

for (nm in names(subsamples)) {
  cat(sprintf("\n--- E: DSEI falsification — %s ---\n", nm))
  sub <- subsamples[[nm]]
  cat(sprintf("  Obs: %d, Polo bases: %d\n", nrow(sub), n_distinct(sub$seqid)))

  iv_cols <- grep(paste0(gp_suffix, "$"), names(sub), value = TRUE)
  x_mat <- as.matrix(sub[, iv_cols])
  f <- sub$seqid
  xr <- lm(x_mat ~ as.factor(f))$residuals
  xr_var <- apply(xr, 2, var)
  threshold <- quantile(xr_var, 0.45)
  selected <- names(xr_var[xr_var > threshold])
  cat(sprintf("  Threshold (q=0.45): %.1f | Passing: %d of %d\n",
              threshold, length(selected), length(xr_var)))

  # PostFE output: only goldmine_area + seqid + year + selected IVs
  df_postfe <- data.frame(
    goldmine_area = sub$goldmine_area,
    seqid = sub$seqid,
    year = sub$year
  )
  df_postfe <- cbind(df_postfe, sub[, selected])
  write_csv(df_postfe, paste0(out_dir, "dsei_falsification_", nm, "_postFE.csv"))
  cat(sprintf("  Saved: %d rows x %d cols\n", nrow(df_postfe), ncol(df_postfe)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Post-FE generation complete ===\n")
postfe_files <- list.files(out_dir, pattern = "_postFE\\.csv$")
cat(sprintf("Generated %d post-FE files:\n", length(postfe_files)))
for (f in sort(postfe_files)) {
  n_iv <- length(grep("2yrgp", readLines(paste0(out_dir, f), n = 1)))
  cat(sprintf("  %-65s %s\n", f, format(file.size(paste0(out_dir, f)), big.mark = ",")))
}
