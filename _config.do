version 18
set more off
set linesize 255

local env_root : environment EVOLVE_HBV_ROOT
if "`env_root'" != "" {
    global PROJ "`env_root'"
}
else {
global PROJ "`c(pwd)'"
}

if "${EVOLVE_HBV_DATA_RAW}" != "" {
    global DATA "$EVOLVE_HBV_DATA_RAW"
}
else {
    local env_data_raw : environment EVOLVE_HBV_DATA_RAW
    if "`env_data_raw'" != "" {
        global DATA "`env_data_raw'"
    }
    else {
        global DATA "$PROJ/data/raw"
        capture confirm file "/Users/lusandamazibuko/Library/CloudStorage/OneDrive-AHRI/Documents/EVOLVE Vukuzazi/data/Vukuzazi_mortality_analysis.dta"
        if !_rc {
            global DATA "/Users/lusandamazibuko/Library/CloudStorage/OneDrive-AHRI/Documents/EVOLVE Vukuzazi/data"
        }
    }
}

global DERIVED "$PROJ/data/derived"
global OUT     "$PROJ/outputs"
global FIGURES "$OUT/figures"
global TABLES  "$OUT/tables"
global LOGS    "$OUT/logs"

foreach d in "$DATA" "$DERIVED" "$OUT" "$FIGURES" "$TABLES" "$LOGS" {
    capture mkdir "`d'"
}

capture program drop drop_stata_tempvars
program define drop_stata_tempvars
    capture ds __0*
    if !_rc & "`r(varlist)'" != "" {
        drop `r(varlist)'
    }
end

capture program drop require_file
program define require_file
    args path
    capture confirm file "`path'"
    if _rc {
        di as error "Required file not found: `path'"
        exit 601
    }
end
