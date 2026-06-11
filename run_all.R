#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Project: Amazon IPLC Mining & Malaria — Reproduction Package
# Purpose: Master script — runs all analysis scripts in order.
#          Reproduces every table, figure, and inline statistic.
#          Sources scripts 02-07 sequentially.
#
# USAGE: From the reproduction/ directory, run:
#   Rscript run_all.R
#
# PREREQUISITES:
#   - R 4.4+ with packages: tidyverse, fixest, modelsummary, kableExtra,
#     gtsummary, janitor, ManyIV, patchwork, sf, ggnewscale
#   - Data files in data/ (included with the package)
#   - Optional: Stata 17+ for LASSO verification scripts in code/iv_selection/
#
# NOTE: The data/ directory is included with the package. No separate
#       data preparation step is needed.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cat("================================================================\n")
cat("  REPRODUCTION PACKAGE: Amazon IPLC Mining & Malaria\n")
cat("================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

t0 <- Sys.time()

# Verify we're in the right directory
if (!file.exists("code/01_setup.R")) {
  stop("Please run this script from the reproduction/ directory.\n",
       "  Expected: cd reproduction && Rscript run_all.R")
}

# Verify data directory exists
if (!dir.exists("data/panels")) {
  stop("Data directory not found. The data/ directory should be included with the package.")
}

# --- Script 02: Municipality analysis ---
cat("\n", strrep("=", 60), "\n")
cat("  02_analysis_municipality.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/02_analysis_municipality.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Script 03: DSEI analysis ---
cat("\n", strrep("=", 60), "\n")
cat("  03_analysis_dsei.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/03_analysis_dsei.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Script 04: Spillover analysis ---
cat("\n", strrep("=", 60), "\n")
cat("  04_analysis_spillover.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/04_analysis_spillover.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Script 05: Extended panel ---
cat("\n", strrep("=", 60), "\n")
cat("  05_analysis_extended_panel.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/05_analysis_extended_panel.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Script 06: IV threshold sensitivity (Figure S5) ---
cat("\n", strrep("=", 60), "\n")
cat("  06_iv_sensitivity.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/06_iv_sensitivity.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Script 07: Combined outputs (Table 1 + all figures) ---
cat("\n", strrep("=", 60), "\n")
cat("  07_combined_outputs.R\n")
cat(strrep("=", 60), "\n")
t1 <- Sys.time()
source("code/07_combined_outputs.R")
cat("  Time:", round(difftime(Sys.time(), t1, units = "mins"), 1), "min\n")

# --- Summary ---
total_time <- round(difftime(Sys.time(), t0, units = "mins"), 1)

cat("\n", strrep("=", 60), "\n")
cat("  REPRODUCTION COMPLETE\n")
cat(strrep("=", 60), "\n")
cat("Total time:", total_time, "min\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# List generated outputs
cat("Tables generated:\n")
tables <- list.files("output/tables/", pattern = "\\.tex$", recursive = TRUE)
for (t in sort(tables)) cat("  ", t, "\n")

cat("\nFigures generated:\n")
figs <- list.files("output/figures/", pattern = "\\.(png|pdf)$", recursive = TRUE)
for (f in sort(figs)) cat("  ", f, "\n")
