version 18
clear all

local thisfile "`c(filename)'"
if "`thisfile'" != "" {
    local suffix "/src/stata/00_run_all.do"
    if substr("`thisfile'", -strlen("`suffix'"), .) == "`suffix'" {
        local root = substr("`thisfile'", 1, strlen("`thisfile'") - strlen("`suffix'"))
        cd "`root'"
    }
}

do "src/stata/_config.do"
do "src/stata/00_preflight.do"

capture log close _all
log using "$LOGS/run_all.smcl", replace

do "src/stata/08_survey_weights.do"
do "src/stata/01_phase1_main_analysis.do"
do "src/stata/02_diagram_numbers.do"
do "src/stata/03_immunity_threshold_sensitivity.do"
do "src/stata/04_threshold_10.do"
do "src/stata/05_immunity_analysis.do"
do "src/stata/06_discordance_analysis.do"
do "src/stata/07_treatment_eligibility.do"

log close
