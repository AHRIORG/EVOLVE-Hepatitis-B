/**********************************************************************
 FIGURE NUMBERS FOR EVOLVE HBV DIAGRAM
 Author : Lusanda Mazibuko
 - Unweighted n and weighted N
 - Exports to: "$OUT/diagram_numbers.xlsx"
**********************************************************************/

version 18
do "src/stata/_config.do"
use "$DERIVED/vk_evolve.dta", clear

/**********************************************************************
 1) Ensure hbvdob exists (0/1/2)
**********************************************************************/
capture confirm variable hbvdob
if _rc {
    capture confirm variable birthyear
    if _rc {
        di as err "ERROR: hbvdob and birthyear both not found."
        exit 198
    }
    gen byte hbvdob = .
    replace hbvdob = 0 if birthyear < 1995
    replace hbvdob = 1 if inrange(birthyear,1995,1999)
    replace hbvdob = 2 if inrange(birthyear,2000,2005)
}
capture label define hbvdob_lbl 0 "Pre-vaccine (<1995)" 1 "Peri-vaccine (1995-1999)" 2 "Post-vaccine (2000-2005)", replace
capture label values hbvdob hbvdob_lbl

/**********************************************************************
 2) Ensure hivstat exists and is labelled
**********************************************************************/
capture confirm variable hivstat
if _rc {
    di as err "ERROR: hivstat not found."
    exit 198
}
capture label define hiv_lbl 0 "HIV-negative" 1 "PLWH", replace
capture label values hivstat hiv_lbl

/**********************************************************************
 3) Clean serology interpretation strings
**********************************************************************/
capture drop _hbsag_i _hbsab_i _ahbc_i
capture drop _miss_hbsag _miss_hbsab _miss_ahbc
capture drop has_hbsag has_hbsab has_ahbc has_all3

gen str40 _hbsag_i = lower(trim(hbsag_interpretation))
gen str40 _hbsab_i = lower(trim(hbsab_interpretation))
gen str40 _ahbc_i  = lower(trim(anti_hbcii_interpretation))

gen byte _miss_hbsag = inlist(_hbsag_i,"","insufficient","equivocal")
gen byte _miss_hbsab = inlist(_hbsab_i,"","insufficient","equivocal")
gen byte _miss_ahbc  = inlist(_ahbc_i ,"","insufficient","equivocal")

gen byte has_hbsag = (_miss_hbsag == 0) if !missing(_hbsag_i)
gen byte has_hbsab = (_miss_hbsab == 0) if !missing(_hbsab_i)
gen byte has_ahbc  = (_miss_ahbc  == 0) if !missing(_ahbc_i)
gen byte has_all3  = (has_hbsag==1 & has_hbsab==1 & has_ahbc==1) ///
                     if !missing(has_hbsag, has_hbsab, has_ahbc)

/**********************************************************************
 4) Create serology outcome categories
**********************************************************************/
capture drop active_disease
gen byte active_disease = .
replace active_disease = 1 if _hbsag_i=="reactive"
replace active_disease = 0 if _hbsag_i=="nonreactive"

capture drop hbvexpo_clear
gen byte hbvexpo_clear = .
replace hbvexpo_clear = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="reactive"
replace hbvexpo_clear = 0 if has_hbsag==1 & has_ahbc==1 & missing(hbvexpo_clear)

capture drop vax_med_immune
gen byte vax_med_immune = .
replace vax_med_immune = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="reactive"
replace vax_med_immune = 0 if has_all3==1 & missing(vax_med_immune)

capture drop susceptible
gen byte susceptible = .
replace susceptible = 1 if _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="nonreactive"
replace susceptible = 0 if has_all3==1 & missing(susceptible)

capture drop sero_cat
gen byte sero_cat = .
replace sero_cat = 1 if _hbsag_i=="reactive"
replace sero_cat = 2 if missing(sero_cat) & _hbsag_i=="nonreactive" & _ahbc_i=="reactive"
replace sero_cat = 3 if missing(sero_cat) & _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="reactive"
replace sero_cat = 4 if missing(sero_cat) & _hbsag_i=="nonreactive" & _ahbc_i=="nonreactive" & _hbsab_i=="nonreactive"

label define serocat_lbl ///
    1 "HBV infection (HBsAg+)" ///
    2 "HBV exposed & cleared (HBsAg-/anti-HBc+)" ///
    3 "Vaccine-mediated immunity (HBsAg-/anti-HBc-/anti-HBs+)" ///
    4 "Susceptible (triple -)", replace
label values sero_cat serocat_lbl

/**********************************************************************
 5) Open postfile to collect diagram numbers
**********************************************************************/
capture postclose DIAG
capture erase "$OUT/diagram_numbers.dta"

postfile DIAG ///
    str40 section ///
    str80 item    ///
    str30 subgroup ///
    long  n       ///
    double n_wt   ///
    using "$OUT/diagram_numbers.dta", replace

/**********************************************************************
 6) Macro to post n and weighted total for a given condition
**********************************************************************/
capture program drop postN
program define postN
    syntax, SECTION(string) ITEM(string) SUBGROUP(string) [IFCOND(string)]

    if "`ifcond'" != "" {
        quietly count if `ifcond'
    }
    else {
        quietly count
    }
    local N = r(N)

    local W = .
    capture confirm variable totwt
    if !_rc {
        if "`ifcond'" != "" {
            quietly summarize totwt if `ifcond'
        }
        else {
            quietly summarize totwt
        }
        local W = r(sum)
    }

    post DIAG ("`section'") ("`item'") ("`subgroup'") (`N') (`W')
end

/**********************************************************************
 7) Sampling frame counts
**********************************************************************/
postN, section("Sampling frame") item("Total enrolled") subgroup("All")
postN, section("Sampling frame") item("Total enrolled") subgroup("HIV-negative") ifcond("hivstat==0")
postN, section("Sampling frame") item("Total enrolled") subgroup("PLWH")         ifcond("hivstat==1")

local vlab : value label hbvdob
foreach c in 0 1 2 {
    local clab "`c'"
    if "`vlab'" != "" {
        local clab : label `vlab' `c'
    }
    postN, section("Sampling frame") item("`clab'") subgroup("Total")        ifcond("hbvdob==`c'")
    postN, section("Sampling frame") item("`clab'") subgroup("HIV-negative") ifcond("hbvdob==`c' & hivstat==0")
    postN, section("Sampling frame") item("`clab'") subgroup("PLWH")         ifcond("hbvdob==`c' & hivstat==1")
}

/**********************************************************************
 8) Serology availability counts
**********************************************************************/
postN, section("Serology availability") item("Has HBsAg result")         subgroup("All") ifcond("has_hbsag==1")
postN, section("Serology availability") item("Has anti-HBs result")       subgroup("All") ifcond("has_hbsab==1")
postN, section("Serology availability") item("Has anti-HBc result")       subgroup("All") ifcond("has_ahbc==1")
postN, section("Serology availability") item("Has all three markers")     subgroup("All") ifcond("has_all3==1")
postN, section("Serology availability") item("Missing one or more markers") subgroup("All") ifcond("has_all3!=1")

/**********************************************************************
 9) Serology outcome counts - overall and by HIV status
**********************************************************************/
foreach hiv in "" "& hivstat==0" "& hivstat==1" {
    if "`hiv'" == ""            local sg "All"
    if "`hiv'" == "& hivstat==0" local sg "HIV-negative"
    if "`hiv'" == "& hivstat==1" local sg "PLWH"

    postN, section("Serology outcomes") item("HBV infection (HBsAg+)")           subgroup("`sg'") ifcond("sero_cat==1 `hiv'")
    postN, section("Serology outcomes") item("HBV exposed and cleared")           subgroup("`sg'") ifcond("sero_cat==2 `hiv'")
    postN, section("Serology outcomes") item("Vaccine-mediated immunity")         subgroup("`sg'") ifcond("sero_cat==3 `hiv'")
    postN, section("Serology outcomes") item("Susceptible to HBV")               subgroup("`sg'") ifcond("sero_cat==4 `hiv'")
    postN, section("Serology outcomes") item("Missing/indeterminate sero category") subgroup("`sg'") ifcond("missing(sero_cat) `hiv'")
}

/**********************************************************************
 10) Close, label and export
**********************************************************************/
postclose DIAG

use "$OUT/diagram_numbers.dta", clear

label var section  "Diagram section"
label var item     "Diagram item"
label var subgroup "Subgroup"
label var n        "Unweighted n"
label var n_wt     "Weighted total (sum of totwt)"

sort section item subgroup
order section item subgroup n n_wt

list, sepby(section) noobs abbreviate(40)

export excel using "$OUT/diagram_numbers.xlsx", firstrow(variables) replace

di as text _newline "Done. Results saved to $OUT/diagram_numbers.xlsx"
