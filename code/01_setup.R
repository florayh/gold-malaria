#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Package loading, path configuration, and shared helper functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Required packages -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(fixest)
library(modelsummary)
library(kableExtra)
library(gtsummary)
library(janitor)
library(ManyIV)
library(patchwork)
library(sf)
library(ggnewscale)

# Check package versions (for reproducibility documentation)
cat("=== Package versions ===\n")
cat("R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
for (pkg in c("tidyverse", "fixest", "modelsummary", "kableExtra",
              "gtsummary", "janitor", "ManyIV", "patchwork", "sf", "ggnewscale")) {
  cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
}
cat("========================\n\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Path configuration -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All paths relative to reproduction/ directory
data_dir    <- "data/panels/"
spill_dir   <- "data/spillover/"
inter_dir   <- "data/intermediate/"
postfe_dir  <- "data/postfe/"
spatial_dir <- "data/spatial/"
table_dir   <- "output/tables/"
figure_dir  <- "output/figures/"

# Create output directories if they don't exist
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Instrument specifications -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Raw rlasso selections (from Stata cluster LASSO, Jan 25, 2026):
#   Municipality (main + weighted): pp3cc_2yrgp pp3_delta_inc_2yrgp pp3_gamma_mlp_2yrgp
#   DSEI (main):    mp3_delta_cs_2yrgp pp3_alfa_bj_2yrgp pp3_delta_in_2yrgp pp3_gamma_mlp_2yrgp c2i_2yrgp pp4_gamma_po_2yrgp
#   DSEI (weighted): mp3_delta_cs_2yrgp pp3_alfa_bj_2yrgp pp3_delta_in_2yrgp pp3_gamma_mlp_2yrgp pp3cc_2yrgp c2i_2yrgp pp4_gamma_po_2yrgp
#
# After backward elimination (drop collinear, then iteratively drop p > 0.10):
#   Municipality: pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp  (pp3cc dropped, p=0.77)
#   DSEI main:    pp3_delta_in_2yrgp + pp3_gamma_mlp_2yrgp   (4 dropped: collinear/insignificant)
#   DSEI weighted: mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp + c2i_2yrgp + pp3_alfa_bj_2yrgp  (3 dropped)
#
# The screening process is run in 02_analysis_municipality.R and 03_analysis_dsei.R
# to verify these selections. See screen_instruments() at end of this file.

# Municipality: Gabro Serra Comprida x GP, Granito Pepita x GP
mu_rlasso_raw <- c("pp3cc_2yrgp", "pp3_delta_inc_2yrgp", "pp3_gamma_mlp_2yrgp")
mu_instruments <- "pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp"

# DSEI: raw rlasso selections and final instruments
dsei_rlasso_raw <- c("mp3_delta_cs_2yrgp", "pp3_alfa_bj_2yrgp", "pp3_delta_in_2yrgp",
                     "pp3_gamma_mlp_2yrgp", "c2i_2yrgp", "pp4_gamma_po_2yrgp")
dsei_instruments <- "pp3_delta_in_2yrgp + pp3_gamma_mlp_2yrgp"

dsei_rlasso_raw_weighted <- c("mp3_delta_cs_2yrgp", "pp3_alfa_bj_2yrgp", "pp3_delta_in_2yrgp",
                              "pp3_gamma_mlp_2yrgp", "pp3cc_2yrgp", "c2i_2yrgp", "pp4_gamma_po_2yrgp")
dsei_weighted_instruments <- "mp3_delta_cs_2yrgp + pp3_delta_in_2yrgp + c2i_2yrgp + pp3_alfa_bj_2yrgp"

# Municipality covariates
mu_select_cov <- "forest_area_km + agriculture_area_km + pop + gdp_per_capita + mean_temp + annual_precip"

# DSEI covariates
dsei_select_cov <- "population_polobase + mean_temp + annual_total_precip + forest"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# fixest dictionary: FE and clustering labels in modelsummary tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setFixest_dict(c(
  cd_mun    = "Mun",
  code_muni = "Mun",
  seqid     = "Subdist.",
  year      = "Year",
  sigla_uf  = "State"
))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Variable labels for modelsummary coef_rename -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coef_rename_mu <- c(
  "fit_gold_mining_area"          = "Gold mining area (km\u00B2)",
  "fit_gold_mining_area_nb100km"  = "Neighbor gold mining area (km\u00B2)",
  "forest_area_km"      = "Forest area (km\u00B2)",
  "agriculture_area_km" = "Agriculture area (km\u00B2)",
  "pop"                 = "Population",
  "gdp_per_capita"      = "GDP per capita",
  "mean_temp"           = "Mean temperature",
  "annual_precip"       = "Annual precipitation"
)

coef_rename_dsei <- c(
  "fit_goldmine_area"    = "Gold mining area (km\u00B2)",
  "population_polobase"  = "Population",
  "mean_temp"            = "Mean temperature",
  "annual_total_precip"  = "Annual precipitation",
  "forest"               = "Forest area (km\u00B2)"
)

# Municipality balance table labels (for gtsummary)
var_labels_mu <- c(
  "pop"                  = "Population",
  "pop_indigenous_2010"  = "Indigenous population (2010 census)",
  "area_km2"             = "Area (km\u00B2)",
  "percent_forest"       = "Forest cover (%)",
  "percent_soy"          = "Soy cover (%)",
  "percent_ag"           = "Agriculture cover (%)",
  "defor_rate"           = "Deforestation rate (%)",
  "health_expenditure"   = "Health expenditure (R$)",
  "health_expenditure_per_capita" = "Health expenditure per capita (R$)"
)

# DSEI balance table labels (for gtsummary)
var_labels_dsei <- c(
  "population_polobase"  = "Population",
  "area_km2"             = "Area (km\u00B2)",
  "mean_temp"            = "Mean temperature",
  "annual_total_precip"  = "Annual precipitation (mm)",
  "percent_forest"       = "Forest cover (%)",
  "percent_soy"          = "Soy cover (%)",
  "percent_all_ag"       = "Agriculture cover (%)",
  "percent_pasture"      = "Pasture cover (%)"
)

# Group labels for balance tables
group_labels_mu <- c(
  "gold_producing" = "Gold-producing",
  "gold_neighbor"  = "Neighbor"
)

group_labels_dsei <- c(
  "gold_producing" = "Gold-producing",
  "gold_neighbor"  = "Neighbor"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# glance_custom.fixest — Extract first-stage F-stat for modelsummary -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This function is automatically called by modelsummary when processing
# fixest objects. It adds the first-stage F-statistic as a GOF row.

glance_custom.fixest <- function(x, ...) {
  out <- data.frame(row.names = 1)

  # Add first-stage F-statistic for IV models
  if (!is.null(x$iv) && x$iv) {
    fs <- fitstat(x, type = "ivf")
    f_stat <- fs[[1]]$stat
    out$`First-stage F` <- sprintf("%.3f", f_stat)
  }

  # Relabel vcov.type value using fixest dict
  dict <- getFixest_dict()
  vcov_fmla <- x$summary_flags$vcov
  if (!is.null(vcov_fmla)) {
    vcov_str <- deparse(vcov_fmla)
    vcov_str <- gsub("~", "", vcov_str)
    vcov_str <- gsub("cluster\\s*", "", vcov_str)
  } else if (!is.null(x$call$cluster)) {
    vcov_str <- as.character(x$call$cluster)
    vcov_str <- gsub("~", "", vcov_str)
    vcov_str <- gsub("cluster\\s*", "", vcov_str)
  } else {
    vcov_str <- NULL
  }
  if (!is.null(vcov_str)) {
    for (nm in names(dict)[order(-nchar(names(dict)))]) {
      vcov_str <- gsub(nm, dict[nm], vcov_str, fixed = TRUE)
    }
    out$vcov.type <- vcov_str
  }

  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standard GOF maps for modelsummary tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gof_map_iv <- tribble(
  ~raw,           ~clean,            ~fmt,
  "nobs",         "Num.Obs.",        0,
  "r.squared",    "R2",              3,
  "First-stage F","First-stage F",   3,
  "FE: cd_mun",   "FE: Mun",        NA,
  "FE: year",     "FE: Year",       NA,
  "vcov.type",    "Std.Errors",      NA
)

gof_map_ols <- tribble(
  ~raw,           ~clean,            ~fmt,
  "nobs",         "Num.Obs.",        0,
  "r.squared",    "R2",              3,
  "FE: cd_mun",   "FE: Mun",        NA,
  "FE: year",     "FE: Year",       NA,
  "vcov.type",    "Std.Errors",      NA
)

gof_map_iv_dsei <- tribble(
  ~raw,           ~clean,            ~fmt,
  "nobs",         "Num.Obs.",        0,
  "r.squared",    "R2",              3,
  "First-stage F","First-stage F",   3,
  "FE: seqid",    "FE: Subdist.", NA,
  "FE: year",     "FE: Year",       NA,
  "vcov.type",    "Std.Errors",      NA
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# fmt_smart — Auto-detect very small values and use more digits -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmt_smart <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || !is.numeric(val)) return(formatC(val, format = "f", digits = 3))
    digits <- 3
    if (abs(val) > 0 & abs(val) < 0.005) digits <- 4
    if (abs(val) > 0 & abs(val) < 0.0005) digits <- 5
    formatC(val, format = "f", digits = digits, big.mark = ",")
  })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# make_iv_summary_table — First-stage diagnostics summary table -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_iv_summary_table <- function(models, spec_names, output_path = NULL) {

  diagnostics <- lapply(seq_along(models), function(i) {
    m <- models[[i]]

    if (is.null(m$iv) || !m$iv) {
      return(data.frame(
        Specification = spec_names[i],
        N = nobs(m),
        `First-stage F` = NA,
        `F p-value` = NA,
        `Sargan p` = NA,
        `N instruments` = NA,
        check.names = FALSE
      ))
    }

    fs <- fitstat(m, type = "ivf")
    f_stat <- round(fs[[1]]$stat, 2)
    f_pval <- fs[[1]]$p

    sargan <- tryCatch({
      s <- fitstat(m, type = "sargan")
      round(s[[1]]$p, 3)
    }, error = function(e) NA)

    n_iv <- length(m$iv_inst_names)

    data.frame(
      Specification = spec_names[i],
      N = nobs(m),
      `First-stage F` = f_stat,
      `F p-value` = format(f_pval, scientific = TRUE, digits = 2),
      `Sargan p` = ifelse(is.na(sargan), "--", as.character(sargan)),
      `N instruments` = n_iv,
      check.names = FALSE
    )
  })

  summary_df <- do.call(rbind, diagnostics)

  if (!is.null(output_path)) {
    kbl <- kableExtra::kbl(summary_df, format = "latex", booktabs = TRUE)
    kableExtra::save_kable(kbl, output_path)
    message("Saved diagnostics table to: ", output_path)
  }

  return(summary_df)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# strip_table_float — Remove table environment wrappers -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strip_table_float <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("strip_table_float: file not found: ", filepath)
    return(invisible(NULL))
  }

  lines <- readLines(filepath, warn = FALSE)

  remove_patterns <- c(
    "^\\s*\\\\begin\\{table\\}",
    "^\\s*\\\\end\\{table\\}",
    "^\\s*\\\\caption\\{",
    "^\\s*\\\\label\\{tab:",
    "^\\s*\\\\centering\\s*$",
    "^\\s*\\\\begin\\{threeparttable\\}",
    "^\\s*\\\\end\\{threeparttable\\}"
  )

  keep <- !grepl(paste(remove_patterns, collapse = "|"), lines)
  new_lines <- lines[keep]

  # Replace minipage wrappers with centering
  new_lines <- gsub("^(\\s*)\\\\begin\\{minipage\\}.*$", "\\1\\\\begin{center}", new_lines)
  new_lines <- gsub("^(\\s*)\\\\end\\{minipage\\}",     "\\1\\\\end{center}",   new_lines)

  # Strip conflicting font size commands
  new_lines <- new_lines[!grepl("^\\s*\\\\fontsize\\{.*\\}\\{.*\\}\\\\selectfont", new_lines)]
  new_lines <- new_lines[!grepl("^\\s*\\\\(footnotesize|small|scriptsize|normalsize)\\s*$", new_lines)]

  # Convert tabular* to tabular
  new_lines <- gsub(
    "\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}(.*?)\\}",
    "\\\\begin{tabular}{\\1}",
    new_lines
  )
  new_lines <- gsub("\\\\end\\{tabular\\*\\}", "\\\\end{tabular}", new_lines)

  # Fix < rendering
  new_lines <- gsub("<0.001", "$<$0.001", new_lines, fixed = TRUE)


  # Wrap in {\small ... }
  new_lines <- c("{\\small", new_lines, "}")

  writeLines(new_lines, filepath)
  message("strip_table_float: processed ", basename(filepath))
  invisible(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# add_spanning_headers — Spanning column headers for tabularray tables -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

add_spanning_headers <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("add_spanning_headers: file not found: ", filepath)
    return(invisible(NULL))
  }

  lines <- readLines(filepath, warn = FALSE)

  inner_close_idx <- grep("tabularray inner close", lines)
  if (length(inner_close_idx) == 0) {
    message("add_spanning_headers: no tabularray inner close found in ", basename(filepath))
    return(invisible(NULL))
  }

  header_idx <- NA
  for (i in (inner_close_idx[1] + 1):length(lines)) {
    if (grepl("&.*\\\\\\\\", lines[i])) {
      header_idx <- i
      break
    }
  }
  if (is.na(header_idx)) {
    message("add_spanning_headers: no header row found in ", basename(filepath))
    return(invisible(NULL))
  }

  header_line <- lines[header_idx]
  header_content <- sub("\\s*\\\\\\\\\\s*$", "", header_line)
  cols <- strsplit(header_content, "\\s*&\\s*")[[1]]
  label_col <- cols[1]
  data_cols <- trimws(cols[-1])

  groups <- list()
  i <- 1
  while (i <= length(data_cols)) {
    name <- data_cols[i]
    count <- 1
    while (i + count <= length(data_cols) && data_cols[i + count] == name) {
      count <- count + 1
    }
    groups <- c(groups, list(list(name = name, count = count)))
    i <- i + count
  }

  has_groups <- any(sapply(groups, function(g) g$count) > 1)
  if (!has_groups) {
    message("add_spanning_headers: no duplicate column headers found in ", basename(filepath))
    return(invisible(NULL))
  }

  n_data_cols <- length(data_cols)
  spanning_parts <- c(label_col)
  numbered_parts <- c(label_col)
  col_num <- 1

  for (g in groups) {
    name <- g$name
    n <- g$count
    display_name <- gsub("Non-falciparum API", "Non-falc. API", name)

    if (n > 1) {
      spanning_parts <- c(spanning_parts, paste0("\\SetCell[c=", n, "]{c} ", display_name))
      spanning_parts <- c(spanning_parts, rep("", n - 1))
    } else {
      spanning_parts <- c(spanning_parts, display_name)
    }

    for (j in 1:n) {
      numbered_parts <- c(numbered_parts, paste0("(", col_num, ")"))
      col_num <- col_num + 1
    }
  }

  spanning_line <- paste0(paste(spanning_parts, collapse = " & "), "  \\\\")
  numbered_line <- paste0(paste(numbered_parts, collapse = " & "), "  \\\\")

  lines <- c(lines[1:(header_idx - 1)],
             spanning_line,
             numbered_line,
             lines[(header_idx + 1):length(lines)])

  for (j in seq_along(lines)) {
    if (grepl("hline\\{\\d+\\}", lines[j])) {
      line <- lines[j]
      matches <- gregexpr("hline\\{(\\d+)\\}", line, perl = TRUE)[[1]]
      if (matches[1] != -1) {
        starts <- as.integer(matches)
        lengths <- attr(matches, "match.length")
        for (k in rev(seq_along(starts))) {
          full_match <- substr(line, starts[k], starts[k] + lengths[k] - 1)
          num <- as.integer(sub("hline\\{(\\d+)\\}", "\\1", full_match))
          if (num >= 2) {
            replacement <- paste0("hline{", num + 1, "}")
            line <- paste0(substr(line, 1, starts[k] - 1),
                          replacement,
                          substr(line, starts[k] + lengths[k], nchar(line)))
          }
        }
        lines[j] <- line
      }
    }
  }

  n_total_cols <- n_data_cols + 1
  inner_open_idx <- grep("tabularray inner open", lines)
  inner_close_idx2 <- grep("tabularray inner close", lines)
  if (length(inner_open_idx) > 0 && length(inner_close_idx2) > 0) {
    hline_new <- paste0("hline{2}={1-", n_total_cols, "}{solid, black, 0.05em},")
    lines <- c(lines[1:(inner_close_idx2[1] - 1)],
               hline_new,
               lines[inner_close_idx2[1]:length(lines)])
  }

  writeLines(lines, filepath)
  message("add_spanning_headers: processed ", basename(filepath))
  invisible(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: build_iv_formula -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build_mu_iv_formula <- function(outcome, covariates = "none") {
  if (covariates == "none") {
    as.formula(paste0(outcome, " ~ 1 | cd_mun + year | gold_mining_area ~ ", mu_instruments))
  } else if (covariates == "select") {
    as.formula(paste0(outcome, " ~ ", mu_select_cov, " | cd_mun + year | gold_mining_area ~ ", mu_instruments))
  }
}

build_dsei_iv_formula <- function(outcome, covariates = "none") {
  if (covariates == "none") {
    as.formula(paste0(outcome, " ~ 1 | seqid + year | goldmine_area ~ ", dsei_instruments))
  } else if (covariates == "select") {
    as.formula(paste0(outcome, " ~ ", dsei_select_cov, " | seqid + year | goldmine_area ~ ", dsei_instruments))
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Instrument screening: backward elimination -----
# Mirrors the manual process in the original analysis scripts:
# 1. Start with all rlasso-selected instruments
# 2. Drop collinear instruments (identified by fixest)
# 3. Iteratively drop the least significant instrument (p > threshold)
# 4. Stop when all remaining instruments are individually significant
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screen_instruments <- function(data, treatment, ivs,
                               fe_vars = "cd_mun + year",
                               cluster_var = "cd_mun",
                               weights_var = NULL,
                               p_threshold = 0.10) {
  current_ivs <- ivs
  for (iter in 1:length(ivs)) {
    if (length(current_ivs) == 0) break
    iv_str <- paste(current_ivs, collapse = " + ")
    fmla <- as.formula(paste0(treatment, " ~ ", iv_str, " | ", fe_vars))

    if (!is.null(weights_var)) {
      w_fmla <- as.formula(paste0("~ ", weights_var))
      m <- tryCatch(feols(fmla, data = data, cluster = cluster_var, weights = w_fmla),
                    error = function(e) NULL)
    } else {
      m <- tryCatch(feols(fmla, data = data, cluster = cluster_var),
                    error = function(e) NULL)
    }
    if (is.null(m)) { cat("    Screening error at iteration", iter, "\n"); break }

    # Step 1: Drop collinear instruments
    collin_ivs <- m$collin.var[m$collin.var %in% current_ivs]
    if (length(collin_ivs) > 0) {
      cat("    Dropping collinear:", paste(collin_ivs, collapse = ", "), "\n")
      current_ivs <- setdiff(current_ivs, collin_ivs)
      next
    }

    # Step 2: Drop least significant if p > threshold
    ct <- summary(m)$coeftable
    pvals <- setNames(ct[, "Pr(>|t|)"], rownames(ct))
    iv_pvals <- pvals[names(pvals) %in% current_ivs]
    worst <- which.max(iv_pvals)

    if (iv_pvals[worst] > p_threshold) {
      cat("    Dropping insignificant (p=", round(iv_pvals[worst], 3), "):",
          names(worst), "\n")
      current_ivs <- setdiff(current_ivs, names(worst))
    } else {
      break
    }
  }
  cat("    Retained", length(current_ivs), "of", length(ivs), "instruments:",
      paste(current_ivs, collapse = ", "), "\n")
  return(current_ivs)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# compute_rotemberg_weights — Pure-R Rotemberg weight decomposition -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Replicates bartik.weight::bw() (https://github.com/jjchern/bartik.weight)
# without the compiled C++ dependency, which fails to build on R >= 4.5
# due to a C++11/C++14 incompatibility in RcppArmadillo.
# The package is GitHub-only (not on CRAN) and appears unmaintained.
# Formulas: Goldsmith-Pinkham, Sorkin, and Swift (2020, QJE).
#
# Arguments:
#   y     — numeric vector (length n), FE-residualized outcome
#   x     — numeric vector (length n), FE-residualized treatment
#   Z_mat — numeric matrix (n x K), FE-residualized instruments
#   G_vec — numeric vector (length K), global shifter (gold price per instrument-year pair)
#
# Returns: list with alpha (K-vector of Rotemberg weights) and beta (K-vector of
#          just-identified IV estimates per instrument-year pair)

compute_rotemberg_weights <- function(y, x, Z_mat, G_vec) {
  # Demean x and y (replicates M_W projection with intercept-only controls)
  xx <- x - mean(x)
  yy <- y - mean(y)
  # Z_mat is already FE-residualized; passed through unchanged (matches bw.cpp)

  Ztx <- as.vector(crossprod(Z_mat, xx))  # K-vector: Z'xx
  Zty <- as.vector(crossprod(Z_mat, yy))  # K-vector: Z'yy

  denom <- sum(G_vec * Ztx)               # scalar: G'Z'xx
  alpha <- (G_vec * Ztx) / denom          # K-vector
  beta  <- Zty / Ztx                      # K-vector

  list(alpha = alpha, beta = beta)
}

message("\n=== Reproduction package setup complete ===\n")
