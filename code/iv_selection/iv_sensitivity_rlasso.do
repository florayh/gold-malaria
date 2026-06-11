* =============================================================================
* Project: Amazon IPLC Mining & Malaria -- Reproduction Package
* Purpose: Figure S5 step 2 (OPTIONAL): Run rlasso across all post-FE quantile
*          threshold files for IV threshold sensitivity analysis.
*          Reads 28 panels from data/intermediate/iv_sensitivity/,
*          writes data/intermediate/rlasso_sweep_selections.csv.
*
*          OPTIONAL -- pre-computed selections are included. Run
*          iv_sensitivity_postfe.R first to generate the input files.
*
* REQUIRES: Stata 17+, lassopack (ssc install lassopack)
* OUTPUT:   data/intermediate/rlasso_sweep_selections.csv
* =============================================================================

* Run from reproduction/ directory:
*   cd reproduction && stata-se -b do code/iv_selection/iv_sensitivity_rlasso.do

local data_dir "data/intermediate/iv_sensitivity"

* =============================================================================
* OUTPUT CSV SETUP
* =============================================================================
local output_csv "data/intermediate/rlasso_sweep_selections.csv"

tempname fh
file open `fh' using "`output_csv'", write replace
file write `fh' "level,quantile,n_candidates,n_selected,instruments" _n

* =============================================================================
* MUNICIPALITY FILES
* =============================================================================

local mu_labels "q00 q05 q10 q15 q20 q25 q30 q35 q40 q45 q50 q55 q60 dropzeros"

foreach qlabel of local mu_labels {
    di ""
    di "============================================================"
    di "MUNICIPALITY: `qlabel'"
    di "============================================================"

    clear
    capture confirm file "`data_dir'/mu_`qlabel'_postFE.csv"
    if _rc != 0 {
        di "  FILE NOT FOUND: mu_`qlabel'_postFE.csv -- skipping"
        continue
    }

    import delimited "`data_dir'/mu_`qlabel'_postFE.csv"

    * Panel setup
    destring cd_mun, replace
    xtset cd_mun year

    * Year dummies
    tab year, gen(y_)
    unab Y : y_*

    * Auto-detect IVs (all columns ending in _2yrgp)
    capture ds *_2yrgp
    if _rc != 0 {
        di "  No IV candidates found -- skipping"
        file write `fh' "mu,`qlabel',0,0," _n
        continue
    }
    local potentialiv `r(varlist)'
    local n_iv : word count `potentialiv'
    di "  Candidates: `n_iv'"

    * Run rlasso
    rlasso gold_mining_area `Y' `potentialiv', fe robust cluster(cd_mun) partial(`Y')

    local sel = e(selected)
    local n_sel : word count `sel'
    di "  Selected: `n_sel'"
    di "  Instruments: `sel'"

    file write `fh' "mu,`qlabel',`n_iv',`n_sel',`sel'" _n
}


* =============================================================================
* DSEI FILES
* =============================================================================

local dsei_labels "q00 q05 q10 q15 q20 q25 q30 q35 q40 q45 q50 q55 q60 dropzeros"

foreach qlabel of local dsei_labels {
    di ""
    di "============================================================"
    di "DSEI: `qlabel'"
    di "============================================================"

    clear
    capture confirm file "`data_dir'/dsei_`qlabel'_postFE.csv"
    if _rc != 0 {
        di "  FILE NOT FOUND: dsei_`qlabel'_postFE.csv -- skipping"
        continue
    }

    import delimited "`data_dir'/dsei_`qlabel'_postFE.csv"

    * Panel setup
    xtset seqid year

    * Year dummies
    tab year, gen(y_)
    unab Y : y_*

    * Auto-detect IVs
    capture ds *_2yrgp
    if _rc != 0 {
        di "  No IV candidates found -- skipping"
        file write `fh' "dsei,`qlabel',0,0," _n
        continue
    }
    local potentialiv `r(varlist)'
    local n_iv : word count `potentialiv'
    di "  Candidates: `n_iv'"

    * Run rlasso
    rlasso goldmine_area `Y' `potentialiv', fe robust cluster(seqid) partial(`Y')

    local sel = e(selected)
    local n_sel : word count `sel'
    di "  Selected: `n_sel'"
    di "  Instruments: `sel'"

    file write `fh' "dsei,`qlabel',`n_iv',`n_sel',`sel'" _n
}


* =============================================================================
file close `fh'
di ""
di "============================================================"
di "All rlasso runs complete."
di "Results saved to: `output_csv'"
di "============================================================"
