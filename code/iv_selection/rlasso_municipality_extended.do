* =============================================================================
* Project: Amazon IPLC Mining & Malaria -- Reproduction Package
* Purpose: Cluster LASSO instrument selection for extended municipality panel
*          (2003-2024, 296 municipalities). Verifies that rlasso selects the
*          instruments used in the extended panel specification (Table S26).
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* EXPECTED OUTPUT:
*   pp3cc_2yrgp mp3_delta_cs_2yrgp pp3_delta_inc_2yrgp
*   (all 3 used directly)
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_municipality_extended.do

local gp_version "2yrgp"

* =============================================================================
* IMPORT DATA
* =============================================================================
clear
import delimited "data/postfe/mu_panel_2yrgp_rockarea_20260316_goldneighbor_postFE.csv"

* =============================================================================
* PANEL SETUP
* =============================================================================
xtset cd_mun year

tab year, gen(y_)
unab Y : y_*

* =============================================================================
* AUTO-DETECT INSTRUMENT VARIABLES
* =============================================================================
ds *_`gp_version'
local potentialiv `r(varlist)'

di "=========================================="
di "MUNICIPALITY -- EXTENDED PANEL (2003-2024)"
di "Number of potential instruments: `: word count `potentialiv''"
di "=========================================="

* =============================================================================
* RUN RLASSO
* =============================================================================
rlasso gold_mining_area `Y' `potentialiv', fe robust cluster(cd_mun) partial(`Y')

di ""
di "SELECTED INSTRUMENTS:"
display e(selected)
di ""
di "=========================================="
di "Expected: pp3cc_2yrgp mp3_delta_cs_2yrgp pp3_delta_inc_2yrgp"
di "=========================================="
