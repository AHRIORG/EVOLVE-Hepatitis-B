version 18
clear all

local thisfile "`c(filename)'"
if "`thisfile'" != "" {
    local suffix "/src/stata/00_preflight.do"
    if substr("`thisfile'", -strlen("`suffix'"), .) == "`suffix'" {
        local root = substr("`thisfile'", 1, strlen("`thisfile'") - strlen("`suffix'"))
        cd "`root'"
    }
}

do "src/stata/_config.do"

di as text "Project root: $PROJ"
di as text "Raw data:     $DATA"
di as text "Derived data: $DERIVED"
di as text "Outputs:      $OUT"

require_file "$DATA/Vukuzazi_mortality_analysis.dta"
require_file "$DATA/VUKUZAZI EVOLVE HBV STUDY RESULTS_13MARCH2024_MOD.xlsx"
require_file "$DATA/EvoLVE AHRI_PLA_for HBV_VL-HBeAg-Alt_RESULTS 18 Oct 2024 1 Mod.xlsx"
require_file "$DATA/AHRI_SER_for Alt_RESULTS 2025-01-25 Mod.xlsx"
require_file "$DATA/list_for_EVOLVE_adult_2004_2005.xlsx"
require_file "$DATA/vukuzazi for sampling.dta"
require_file "$DATA/Vukuzazi_weights_2022_age_sex.dta"
require_file "$DATA/Mapping EVOLVE.dta"

capture confirm file "$DATA/EVOLVE - April Data.dta"
if _rc {
    require_file "$DATA/EVOLVE - April Data .dta"
}

di as result "Stata preflight passed."
