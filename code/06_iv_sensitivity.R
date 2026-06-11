#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Figure S5: Read rlasso selections from quantile sweep,
#          run backward elimination + 2SLS at each threshold, save
#          coefficient stability results, and plot Figure S5.
#
#          Reads: data/intermediate/rlasso_sweep_selections.csv (pre-computed)
#          Writes: data/intermediate/iv_sensitivity_results.csv
#                  output/figures/fig_iv_threshold_sensitivity.pdf
#
# Author: Flora He (assisted by Claude Code)
# Start date: 2026-06-11
# Code review: NOT REVIEWED
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("code/01_setup.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Configuration -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rlasso_file <- paste0(inter_dir, "rlasso_sweep_selections.csv")
results_file <- paste0(inter_dir, "iv_sensitivity_results.csv")

# Main-spec chosen quantiles (for marking on the plot)
mu_main_q   <- 0.25
dsei_main_q <- 0.45

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_df <- read_csv(paste0(data_dir, "mu_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
mu_df$cd_mun <- as.character(mu_df$cd_mun)

dsei_df <- read_csv(paste0(data_dir, "dsei_panel_goldneighbor_0319.csv"), show_col_types = FALSE)
dsei_df$seqid <- as.character(dsei_df$seqid)

# Compute parasite-specific API for DSEI
dsei_df$falciparum_api <- 1000 * dsei_df$malaria_n_falciparum_mix / dsei_df$population_polobase
dsei_df$non_f_api <- 1000 * dsei_df$malaria_n_not_f / dsei_df$population_polobase

rlasso_sel <- read_csv(rlasso_file, show_col_types = FALSE)
cat("Loaded rlasso selections:", nrow(rlasso_sel), "rows\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run estimation for each threshold -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results <- tibble()

for (i in seq_len(nrow(rlasso_sel))) {
  row <- rlasso_sel[i, ]
  level <- row$level
  q_label <- row$quantile
  n_cand <- row$n_candidates
  n_sel  <- row$n_selected

  cat(sprintf("\n--- %s %s (%d candidates, %d rlasso-selected) ---\n",
              toupper(level), q_label, n_cand, n_sel))

  # Parse instruments (space-separated in the CSV)
  iv_string <- row$instruments
  if (is.na(iv_string) || iv_string == "" || n_sel == 0) {
    cat("  No instruments selected — recording NA\n")
    results <- bind_rows(results, tibble(
      level = level, q_label = q_label, n_candidates = n_cand,
      n_rlasso = 0, n_after_screen = 0, instruments = NA_character_,
      coef = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      fstat = NA_real_, pval = NA_real_
    ))
    next
  }

  ivs <- str_split(iv_string, "\\s+")[[1]]
  ivs <- ivs[ivs != ""]

  # Select data and parameters based on level
  if (level == "mu") {
    dat <- mu_df
    treatment <- "gold_mining_area"
    outcome <- "malaria_allpop_api"
    unit_var <- "cd_mun"
    select_cov <- mu_select_cov
  } else {
    dat <- dsei_df
    treatment <- "goldmine_area"
    outcome <- "malaria_api_by_polobase"
    unit_var <- "seqid"
    select_cov <- dsei_select_cov
  }

  # Backward elimination (reuses screen_instruments from 01_setup.R)
  cat("  Screening instruments...\n")
  screened_ivs <- screen_instruments(dat, treatment, ivs,
                                     fe_vars = paste0(unit_var, " + year"),
                                     cluster_var = unit_var)

  if (length(screened_ivs) == 0) {
    cat("  All instruments eliminated — recording NA\n")
    results <- bind_rows(results, tibble(
      level = level, q_label = q_label, n_candidates = n_cand,
      n_rlasso = length(ivs), n_after_screen = 0, instruments = NA_character_,
      coef = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      fstat = NA_real_, pval = NA_real_
    ))
    next
  }

  # Run 2SLS with selected covariates
  iv_str <- paste(screened_ivs, collapse = " + ")
  fmla_2sls <- as.formula(paste0(
    outcome, " ~ ", select_cov, " | ", unit_var, " + year | ",
    treatment, " ~ ", iv_str
  ))

  mod <- tryCatch(
    feols(fmla_2sls, data = dat, cluster = unit_var),
    error = function(e) { cat("  2SLS error:", e$message, "\n"); NULL }
  )

  if (is.null(mod)) {
    results <- bind_rows(results, tibble(
      level = level, q_label = q_label, n_candidates = n_cand,
      n_rlasso = length(ivs), n_after_screen = length(screened_ivs),
      instruments = iv_str,
      coef = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      fstat = NA_real_, pval = NA_real_
    ))
    next
  }

  # Extract results
  ct <- summary(mod)$coeftable
  treat_row <- which(rownames(ct) == paste0("fit_", treatment))
  if (length(treat_row) == 0) treat_row <- 1

  beta <- ct[treat_row, "Estimate"]
  se_val <- ct[treat_row, "Std. Error"]
  pv <- ct[treat_row, "Pr(>|t|)"]
  ci_lo <- beta - 1.96 * se_val
  ci_hi <- beta + 1.96 * se_val

  fs <- tryCatch(fitstat(mod, type = "ivf")[[1]]$stat, error = function(e) NA_real_)

  cat(sprintf("  coef=%.3f se=%.3f F=%.1f IVs=%d\n",
              beta, se_val, fs, length(screened_ivs)))

  results <- bind_rows(results, tibble(
    level = level, q_label = q_label, n_candidates = n_cand,
    n_rlasso = length(ivs), n_after_screen = length(screened_ivs),
    instruments = iv_str,
    coef = beta, se = se_val, ci_lo = ci_lo, ci_hi = ci_hi,
    fstat = fs, pval = pv
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Save results -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write_csv(results, results_file)
cat(sprintf("\nResults saved to: %s (%d rows)\n", results_file, nrow(results)))
print(results %>% select(level, q_label, n_candidates, n_rlasso, n_after_screen, coef, se, fstat))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Figure S5: IV threshold sensitivity plot -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== Figure S5: IV threshold sensitivity ===\n")

parse_q <- function(q_label) {
  as.numeric(sub("^q", "", q_label)) / 100
}

make_coef_plot <- function(res, level_name, main_q, title_label, show_xlab = TRUE) {
  df_plot <- res %>%
    filter(level == level_name, !is.na(coef), q_label != "dropzeros") %>%
    mutate(q_num = parse_q(q_label),
           is_main = abs(q_num - main_q) < 0.001)

  if (nrow(df_plot) == 0) return(NULL)

  ggplot(df_plot, aes(x = q_num, y = coef)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = main_q, linetype = "dashed", color = "firebrick3", linewidth = 0.5) +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                    color = "steelblue4", size = 0.4, linewidth = 0.6) +
    geom_point(data = df_plot %>% filter(is_main),
               aes(x = q_num, y = coef),
               color = "firebrick3", size = 3, shape = 18) +
    labs(x = if (show_xlab) "Post-FE variance quantile threshold" else NULL,
         y = "2SLS coefficient",
         title = paste0(title_label, ": coefficient stability")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          panel.grid.minor = element_blank())
}

p_coef_mu <- make_coef_plot(results, "mu", mu_main_q, "Municipality")
p_coef_dsei <- make_coef_plot(results, "dsei", dsei_main_q, "DSEI polo base")

if (!is.null(p_coef_mu) && !is.null(p_coef_dsei)) {
  p_s5 <- p_coef_mu | p_coef_dsei
  ggsave(paste0(figure_dir, "fig_iv_threshold_sensitivity.pdf"),
         p_s5, width = 12, height = 5)
  cat("Saved:", paste0(figure_dir, "fig_iv_threshold_sensitivity.pdf"), "\n")
}
