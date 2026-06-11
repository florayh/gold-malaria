* =============================================================================
* Project: Amazon IPLC Mining & Malaria — Reproduction Package
* Purpose: Cluster LASSO instrument selection for municipality-level analysis
*          (Belloni et al. 2016). Verifies that rlasso selects the instruments
*          used in the main specification.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* EXPECTED OUTPUT:
*   Main (unweighted): pp3cc_2yrgp pp3_delta_inc_2yrgp pp3_gamma_mlp_2yrgp
*     -> pp3cc dropped in R backward elimination (p=0.77)
*     -> Final: pp3_delta_inc_2yrgp + pp3_gamma_mlp_2yrgp
*
*   Weighted (Table S17): pp3_delta_inc_2yrgp pp3_gamma_mlp_2yrgp
*
*   Drop Itaituba+Jacareacanga (Table S13):
*     pp3_delta_inc_2yrgp pp3_gamma_mlp_2yrgp
* =============================================================================

local gp_version "2yrgp"

* =============================================================================
* 1. MAIN SPECIFICATION (unweighted)
* =============================================================================
clear
* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_municipality.do
import delimited "data/postfe/mu_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_postFE.csv"

xtset cd_mun year

tab year, gen(y_)
unab Y : y_*

ds *_`gp_version'
local potentialiv `r(varlist)'

di "=========================================="
di "MUNICIPALITY — UNWEIGHTED (MAIN)"
di "Number of potential instruments: `: word count `potentialiv''"
di "=========================================="

rlasso gold_mining_area `Y' `potentialiv', fe robust cluster(cd_mun) partial(`Y')

di ""
di "SELECTED INSTRUMENTS:"
display e(selected)
di ""

* =============================================================================
* 2. WEIGHTED SPECIFICATION (Table S17)
* =============================================================================
drop y_*

tab year, gen(y_)
unab Y : y_*

di "=========================================="
di "MUNICIPALITY — POPULATION WEIGHTED (Table S17)"
di "=========================================="

rlasso gold_mining_area `Y' `potentialiv' [pweight=pop], fe robust cluster(cd_mun) partial(`Y')

di ""
di "SELECTED INSTRUMENTS (WEIGHTED):"
display e(selected)
di ""

* =============================================================================
* 3. DROP ITAITUBA + JACAREACANGA (Table S13)
* =============================================================================
clear
import delimited "data/postfe/mu_panel_2yrgp_rockarea_20260316_03-19_goldneighbor_dropIJ_postFE.csv"

* Verify exclusion
count if cd_mun == 150360
assert r(N) == 0
count if cd_mun == 150375
assert r(N) == 0

xtset cd_mun year

tab year, gen(y_)
unab Y : y_*

ds *_`gp_version'
local potentialiv `r(varlist)'

di "=========================================="
di "MUNICIPALITY — DROP ITAITUBA + JACAREACANGA (Table S13)"
di "Number of potential instruments: `: word count `potentialiv''"
di "=========================================="

rlasso gold_mining_area `Y' `potentialiv', fe robust cluster(cd_mun) partial(`Y')

di ""
di "SELECTED INSTRUMENTS (DROP IJ):"
display e(selected)
di ""
di "=========================================="
di "VERIFICATION COMPLETE"
di "=========================================="
