#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Figure S5 step 1 (OPTIONAL): Generate post-FE filtered IV files
#          across a grid of quantile thresholds (0.00–0.60) for the IV
#          threshold sensitivity analysis.
#          Output: 28 postFE panels in data/intermediate/iv_sensitivity/
#
#          OPTIONAL — the pre-computed rlasso selections in
#          data/intermediate/rlasso_sweep_selections.csv make this unnecessary.
#          Run this + iv_sensitivity_rlasso.do to regenerate from scratch.
#
# Author: Flora He (assisted by Claude Code)
# Start date: 2026-06-11
# Code review: NOT REVIEWED
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Configuration -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

panels_dir <- "data/panels/"
out_dir    <- "data/intermediate/iv_sensitivity/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gp_suffix <- "_2yrgp"

# Quantile grid: 0.00 to 0.60 in steps of 0.05, plus a "drop zeros only" variant
quantiles <- seq(0, 0.60, by = 0.05)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: post-FE filter at a given quantile -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generate_postfe <- function(df, interaction_cols, endogenous, unit_var,
                            outcome_var, extra_vars, quantile_val,
                            drop_zeros_only = FALSE) {
  x_mat <- as.matrix(df[, interaction_cols])
  f <- df[[unit_var]]
  xr <- lm(x_mat ~ as.factor(f))$residuals
  xr_var <- apply(xr, 2, var)

  if (drop_zeros_only) {
    selected_names <- names(xr_var[xr_var > 0])
  } else if (quantile_val == 0) {
    selected_names <- names(xr_var)
  } else {
    threshold <- quantile(xr_var, quantile_val)
    selected_names <- names(xr_var[xr_var > threshold])
  }

  out <- df %>%
    select(all_of(c(unit_var, "year", endogenous, outcome_var, extra_vars, selected_names)))

  list(data = out, n_ivs = length(selected_names))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Municipality sweep -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message("=== Municipality panel ===")
mu_df <- read_csv(paste0(panels_dir, "mu_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
mu_df$cd_mun <- as.character(mu_df$cd_mun)

mu_interaction_cols <- names(mu_df)[grep(paste0(gp_suffix, "$"), names(mu_df))]
message(paste("Total interaction columns:", length(mu_interaction_cols)))

mu_extra <- c("pop", "forest_area_km", "agriculture_area_km",
              "gdp_per_capita", "mean_temp", "annual_precip")

mu_summary <- tibble(level = character(), q = character(), n_ivs = integer())

for (q in quantiles) {
  q_label <- sprintf("q%02d", round(q * 100))
  message(paste0("  MU ", q_label, " (q=", q, ")..."))
  result <- generate_postfe(mu_df, mu_interaction_cols, "gold_mining_area",
                            "cd_mun", "malaria_allpop_api", mu_extra, q)
  outfile <- file.path(out_dir, paste0("mu_", q_label, "_postFE.csv"))
  write_csv(result$data, outfile)
  mu_summary <- bind_rows(mu_summary, tibble(level = "mu", q = q_label, n_ivs = result$n_ivs))
  message(paste0("    -> ", result$n_ivs, " IVs"))
}

# Drop-zeros-only variant
message("  MU dropzeros...")
result_dz <- generate_postfe(mu_df, mu_interaction_cols, "gold_mining_area",
                             "cd_mun", "malaria_allpop_api", mu_extra, 0,
                             drop_zeros_only = TRUE)
write_csv(result_dz$data, file.path(out_dir, "mu_dropzeros_postFE.csv"))
mu_summary <- bind_rows(mu_summary, tibble(level = "mu", q = "dropzeros", n_ivs = result_dz$n_ivs))
message(paste0("    -> ", result_dz$n_ivs, " IVs"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DSEI sweep -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message("\n=== DSEI panel ===")
dsei_df <- read_csv(paste0(panels_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
dsei_df$seqid <- as.character(dsei_df$seqid)

dsei_interaction_cols <- names(dsei_df)[grep(paste0(gp_suffix, "$"), names(dsei_df))]
message(paste("Total interaction columns:", length(dsei_interaction_cols)))

dsei_extra <- c("population_polobase", "mean_temp", "annual_total_precip", "forest")

dsei_summary <- tibble(level = character(), q = character(), n_ivs = integer())

for (q in quantiles) {
  q_label <- sprintf("q%02d", round(q * 100))
  message(paste0("  DSEI ", q_label, " (q=", q, ")..."))
  result <- generate_postfe(dsei_df, dsei_interaction_cols, "goldmine_area",
                            "seqid", "malaria_api_by_polobase", dsei_extra, q)
  outfile <- file.path(out_dir, paste0("dsei_", q_label, "_postFE.csv"))
  write_csv(result$data, outfile)
  dsei_summary <- bind_rows(dsei_summary, tibble(level = "dsei", q = q_label, n_ivs = result$n_ivs))
  message(paste0("    -> ", result$n_ivs, " IVs"))
}

# Drop-zeros-only variant
message("  DSEI dropzeros...")
result_dz <- generate_postfe(dsei_df, dsei_interaction_cols, "goldmine_area",
                             "seqid", "malaria_api_by_polobase", dsei_extra, 0,
                             drop_zeros_only = TRUE)
write_csv(result_dz$data, file.path(out_dir, "dsei_dropzeros_postFE.csv"))
dsei_summary <- bind_rows(dsei_summary, tibble(level = "dsei", q = "dropzeros", n_ivs = result_dz$n_ivs))
message(paste0("    -> ", result_dz$n_ivs, " IVs"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save summary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_all <- bind_rows(mu_summary, dsei_summary)
write_csv(summary_all, file.path(out_dir, "postfe_sweep_summary.csv"))
message(paste("\nSummary saved. Generated", nrow(summary_all), "postFE files."))
print(summary_all)
