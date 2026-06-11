* =============================================================================
* Project: Amazon IPLC Mining & Malaria — Reproduction Package
* Purpose: Cluster LASSO instrument selection for DSEI polo base analysis
*          (Belloni et al. 2016). Verifies that rlasso selects the instruments
*          used in the main specification.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* EXPECTED OUTPUT:
*   Main (unweighted): pp3_delta_in_2yrgp pp3_gamma_mlp_2yrgp
*     (+ 3 others that are collinear or have lower F; dropped in R)
*
*   Weighted (Table S18):
*     mp3_delta_cs_2yrgp pp3_delta_in_2yrgp c2i_2yrgp pp3_alfa_bj_2yrgp
* =============================================================================

local gp_version "2yrgp"

* =============================================================================
* 1. MAIN SPECIFICATION (unweighted)
* =============================================================================
clear
* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_dsei.do
import delimited "data/postfe/dsei_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_postFE.csv"

xtset seqid year

tab year, gen(y_)
unab Y : y_*

ds *_`gp_version'
local potentialiv `r(varlist)'

di "=========================================="
di "DSEI — UNWEIGHTED (MAIN)"
di "Number of potential instruments: `: word count `potentialiv''"
di "=========================================="

rlasso goldmine_area `Y' `potentialiv', fe robust cluster(seqid) partial(`Y')

di ""
di "SELECTED INSTRUMENTS:"
display e(selected)
di ""

* =============================================================================
* 2. WEIGHTED SPECIFICATION (Table S18)
* =============================================================================
drop y_*

tab year, gen(y_)
unab Y : y_*

di "=========================================="
di "DSEI — POPULATION WEIGHTED (Table S18)"
di "=========================================="

rlasso goldmine_area `Y' `potentialiv' [pweight=population_polobase], fe robust cluster(seqid) partial(`Y')

di ""
di "SELECTED INSTRUMENTS (WEIGHTED):"
display e(selected)
di ""
di "=========================================="
di "VERIFICATION COMPLETE"
di "=========================================="
