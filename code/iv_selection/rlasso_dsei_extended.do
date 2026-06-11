* =============================================================================
* Project: Amazon IPLC Mining & Malaria -- Reproduction Package
* Purpose: Cluster LASSO instrument selection for extended DSEI panel
*          (2003-2024, 106 polo bases). Verifies that rlasso selects the
*          instruments used in the extended panel specification (Table S27).
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
*
* EXPECTED OUTPUT:
*   pp3cc_2yrgp c2i_2yrgp pp3_alfa_ar_2yrgp pp4_gamma_po_2yrgp
*   (all 4 used directly)
* =============================================================================

* Run from reproduction/ directory: cd reproduction && stata-se -b do code/iv_selection/rlasso_dsei_extended.do

local gp_version "2yrgp"

* =============================================================================
* IMPORT DATA
* =============================================================================
clear
import delimited "data/postfe/dsei_panel_2yrgp_rockarea_20260316_goldneighbor_postFE.csv"

* =============================================================================
* PANEL SETUP
* =============================================================================
xtset seqid year

tab year, gen(y_)
unab Y : y_*

* =============================================================================
* AUTO-DETECT INSTRUMENT VARIABLES
* =============================================================================
ds *_`gp_version'
local potentialiv `r(varlist)'

di "=========================================="
di "DSEI -- EXTENDED PANEL (2003-2024)"
di "Number of potential instruments: `: word count `potentialiv''"
di "=========================================="

* =============================================================================
* RUN RLASSO
* =============================================================================
rlasso goldmine_area `Y' `potentialiv', fe robust cluster(seqid) partial(`Y')

di ""
di "SELECTED INSTRUMENTS:"
display e(selected)
di ""
di "=========================================="
di "Expected: pp3cc_2yrgp c2i_2yrgp pp3_alfa_ar_2yrgp pp4_gamma_po_2yrgp"
di "=========================================="
