/*===============================================================*
 | DISCORDANCE ANALYSIS
 | HBsAg-positive / anti-HBc-negative phenotype
 | Purpose:
 |   Reviewer query: absence of anti-HBc among HBsAg-positive cases
 |   Assess whether discordance differs by HBsAg S/CO threshold,
 |   HIV status, HBsAg signal strength, and HBV DNA where available.
 
 
 =================================================================
 | HBsAg-positive / anti-HBc-negative phenotype
 | Includes:
 |   1) Unweighted phenotype counts
 |   2) Weighted population-level discordance estimates
 |   3) Weighted discordance among HBsAg-positive participants
 |   4) Weighted stratified estimates by HIV, sex, birth cohort
 |   5) Clinical/descriptive characterisation: HBsAg S/CO and HBV DNA
 *===============================================================*/

version 18
do "src/stata/_config.do"
use "$DERIVED/vk_evolve.dta", clear
capture drop __0*

/* IMPORTANT:
   Do NOT run: svyset, singleunit(centered)
   That removes the sampling weights.
*/
svyset _n [pweight=totwt], singleunit(centered)
svyset

capture mkdir "$OUT/discordance"

/*---------------------------------------------------------------*
 | 1. Ensure HBsAg numeric S/CO exists
 *---------------------------------------------------------------*/

capture confirm numeric variable hbv_num
if _rc {
    capture confirm numeric variable hbsag_result
    if _rc {
        gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9]+(\.[0-9]+)?)")
    }
    else {
        gen double hbv_num = hbsag_result
    }
}
label var hbv_num "HBsAg S/CO numeric value"

/*---------------------------------------------------------------*
 | 2. Clean anti-HBc interpretation
 *---------------------------------------------------------------*/

cap drop antihbc_clean 

gen str30 antihbc_clean = lower(trim(anti_hbcii_interpretation))

cap drop antihbc_pos
gen byte antihbc_pos = .
replace antihbc_pos = 1 if antihbc_clean=="reactive"
replace antihbc_pos = 0 if antihbc_clean=="nonreactive"

label define antihbc_lbl 0 "anti-HBc negative" 1 "anti-HBc positive", replace
label values antihbc_pos antihbc_lbl
label var antihbc_pos "anti-HBc status"

/*---------------------------------------------------------------*
 | 3. Define HBsAg positivity at both thresholds
 *---------------------------------------------------------------*/

cap drop hbsag_gt1 hbsag_gt10

gen byte hbsag_gt1 = .
replace hbsag_gt1 = 1 if hbv_num > 1  & hbv_num < .
replace hbsag_gt1 = 0 if hbv_num <= 1 & hbv_num < .

gen byte hbsag_gt10 = .
replace hbsag_gt10 = 1 if hbv_num > 10  & hbv_num < .
replace hbsag_gt10 = 0 if hbv_num <= 10 & hbv_num < .

label define hbsag_lbl 0 "HBsAg negative" 1 "HBsAg positive", replace
label values hbsag_gt1 hbsag_lbl
label values hbsag_gt10 hbsag_lbl

label var hbsag_gt1  "HBsAg S/CO >1"
label var hbsag_gt10 "HBsAg S/CO >10"

/*---------------------------------------------------------------*
 | 4. Define discordance among HBsAg-positive individuals
 *---------------------------------------------------------------*/

cap drop discord_gt1 discord_gt10

gen byte discord_gt1 = .
replace discord_gt1 = 1 if hbsag_gt1==1  & antihbc_pos==0
replace discord_gt1 = 0 if hbsag_gt1==1  & antihbc_pos==1

gen byte discord_gt10 = .
replace discord_gt10 = 1 if hbsag_gt10==1 & antihbc_pos==0
replace discord_gt10 = 0 if hbsag_gt10==1 & antihbc_pos==1

label define discord_lbl 0 "HBsAg+/anti-HBc+" 1 "HBsAg+/anti-HBc-", replace
label values discord_gt1 discord_lbl
label values discord_gt10 discord_lbl

label var discord_gt1  "Discordant phenotype among HBsAg S/CO >1"
label var discord_gt10 "Discordant phenotype among HBsAg S/CO >10"

/*---------------------------------------------------------------*
 | 5. Unweighted counts: observed phenotype among HBsAg-positive
 *---------------------------------------------------------------*/

display as text "======================================================"
display as text "UNWEIGHTED DISCORDANCE AMONG HBsAg S/CO >1"
display as text "======================================================"
tab discord_gt1 if hbsag_gt1==1, missing

display as text "======================================================"
display as text "UNWEIGHTED DISCORDANCE AMONG HBsAg S/CO >10"
display as text "======================================================"
tab discord_gt10 if hbsag_gt10==1, missing

/*---------------------------------------------------------------*
 | 6. WEIGHTED population-level prevalence of discordant phenotype
 |    This estimates prevalence among all screened participants.
 *---------------------------------------------------------------*/

cap drop disc_pop_gt1 disc_pop_gt10

gen byte disc_pop_gt1 = .
replace disc_pop_gt1 = 1 if hbsag_gt1==1 & antihbc_pos==0
replace disc_pop_gt1 = 0 if hbsag_gt1==0 | (hbsag_gt1==1 & antihbc_pos==1)

gen byte disc_pop_gt10 = .
replace disc_pop_gt10 = 1 if hbsag_gt10==1 & antihbc_pos==0
replace disc_pop_gt10 = 0 if hbsag_gt10==0 | (hbsag_gt10==1 & antihbc_pos==1)

label var disc_pop_gt1  "Population prevalence: HBsAg S/CO >1 and anti-HBc negative"
label var disc_pop_gt10 "Population prevalence: HBsAg S/CO >10 and anti-HBc negative"

display as text "======================================================"
display as text "WEIGHTED POPULATION PREVALENCE OF DISCORDANT PHENOTYPE"
display as text "======================================================"
svy: mean disc_pop_gt1
svy: mean disc_pop_gt10

/*---------------------------------------------------------------*
 | 7. WEIGHTED proportion discordant among HBsAg-positive
 *---------------------------------------------------------------*/

cap drop sub_gt1 sub_gt10

gen byte sub_gt1  = hbsag_gt1==1  & !missing(discord_gt1)
gen byte sub_gt10 = hbsag_gt10==1 & !missing(discord_gt10)

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE AMONG HBsAg S/CO >1"
display as text "======================================================"
svy, subpop(sub_gt1): mean discord_gt1

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE AMONG HBsAg S/CO >10"
display as text "======================================================"
svy, subpop(sub_gt10): mean discord_gt10

/*---------------------------------------------------------------*
 | 8. WEIGHTED discordance by HIV status
 *---------------------------------------------------------------*/

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY HIV STATUS: S/CO >1"
display as text "======================================================"
tab hivstat discord_gt1 if hbsag_gt1==1, row col missing
svy, subpop(sub_gt1): proportion discord_gt1, over(hivstat)
svy, subpop(sub_gt1): logistic discord_gt1 i.hivstat

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY HIV STATUS: S/CO >10"
display as text "======================================================"
tab hivstat discord_gt10 if hbsag_gt10==1, row col missing
svy, subpop(sub_gt10): proportion discord_gt10, over(hivstat)
svy, subpop(sub_gt10): logistic discord_gt10 i.hivstat

/*---------------------------------------------------------------*
 | 9. WEIGHTED discordance by sex
 *---------------------------------------------------------------*/

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY SEX: S/CO >1"
display as text "======================================================"
tab sex discord_gt1 if hbsag_gt1==1, row col missing
svy, subpop(sub_gt1): proportion discord_gt1, over(sex)
svy, subpop(sub_gt1): logistic discord_gt1 i.sex

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY SEX: S/CO >10"
display as text "======================================================"
tab sex discord_gt10 if hbsag_gt10==1, row col missing
svy, subpop(sub_gt10): proportion discord_gt10, over(sex)
svy, subpop(sub_gt10): logistic discord_gt10 i.sex

/*---------------------------------------------------------------*
 | 10. WEIGHTED discordance by birth cohort
 *---------------------------------------------------------------*/

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY BIRTH COHORT: S/CO >1"
display as text "======================================================"
tab hbvdob discord_gt1 if hbsag_gt1==1, row col missing
svy, subpop(sub_gt1): proportion discord_gt1, over(hbvdob)
svy, subpop(sub_gt1): logistic discord_gt1 i.hbvdob

display as text "======================================================"
display as text "WEIGHTED DISCORDANCE BY BIRTH COHORT: S/CO >10"
display as text "======================================================"
tab hbvdob discord_gt10 if hbsag_gt10==1, row col missing
svy, subpop(sub_gt10): proportion discord_gt10, over(hbvdob)
svy, subpop(sub_gt10): logistic discord_gt10 i.hbvdob

/*---------------------------------------------------------------*
 | 11. HBsAg S/CO signal strength by concordance group
 |     Descriptive clinical characterisation.
 *---------------------------------------------------------------*/

display as text "======================================================"
display as text "HBsAg S/CO DISTRIBUTION BY DISCORDANCE: S/CO >1"
display as text "======================================================"
tabstat hbv_num if hbsag_gt1==1, by(discord_gt1) ///
    stat(n mean sd median iqr min max) columns(statistics)

ranksum hbv_num if hbsag_gt1==1, by(discord_gt1)

display as text "======================================================"
display as text "HBsAg S/CO DISTRIBUTION BY DISCORDANCE: S/CO >10"
display as text "======================================================"
tabstat hbv_num if hbsag_gt10==1, by(discord_gt10) ///
    stat(n mean sd median iqr min max) columns(statistics)

ranksum hbv_num if hbsag_gt10==1, by(discord_gt10)

/*---------------------------------------------------------------*
 | 12. Anti-HBs status among HBsAg-positive participants
 *---------------------------------------------------------------*/

cap drop antihbs_clean antihbs_pos

gen str30 antihbs_clean = lower(trim(hbsab_interpretation))

cap drop antihbs_pos
gen byte antihbs_pos = .
replace antihbs_pos = 1 if antihbs_clean=="reactive"
replace antihbs_pos = 0 if antihbs_clean=="nonreactive"

label define antihbs_lbl 0 "anti-HBs negative" 1 "anti-HBs positive", replace
label values antihbs_pos antihbs_lbl
label var antihbs_pos "anti-HBs status"

display as text "======================================================"
display as text "ANTI-HBs STATUS BY DISCORDANCE: S/CO >1"
display as text "======================================================"
tab antihbs_pos discord_gt1 if hbsag_gt1==1, row col missing
svy, subpop(sub_gt1): proportion discord_gt1, over(antihbs_pos)
svy, subpop(sub_gt1): logistic discord_gt1 i.antihbs_pos

display as text "======================================================"
display as text "ANTI-HBs STATUS BY DISCORDANCE: S/CO >10"
display as text "======================================================"
tab antihbs_pos discord_gt10 if hbsag_gt10==1, row col missing
svy, subpop(sub_gt10): proportion discord_gt10, over(antihbs_pos)
svy, subpop(sub_gt10): logistic discord_gt10 i.antihbs_pos

/*---------------------------------------------------------------*
 | 13. Export individual-level discordance dataset
 *---------------------------------------------------------------*/

preserve

keep individualid iintid hbv_num hbsag_gt1 hbsag_gt10 ///
     antihbc_pos antihbs_pos discord_gt1 discord_gt10 ///
     disc_pop_gt1 disc_pop_gt10 ///
     hivstat sex hbvdob birthyear educ_recoded sescat bmicat ///
     hypertension diabetic onart vl_suppressed totwt

keep if hbsag_gt1==1 | hbsag_gt10==1

capture drop __0*
export excel using "$OUT/discordance/discordance_individual_level.xlsx", ///
    firstrow(variables) replace

capture drop __0*
save "$OUT/discordance/discordance_individual_level.dta", replace

restore

/*---------------------------------------------------------------*
 | 14. Export unweighted summary counts
 *---------------------------------------------------------------*/

preserve

tempfile discsum
postutil clear

postfile D str12 cutoff str30 phenotype double n pct using `discsum', replace

foreach cut in gt1 gt10 {
    quietly count if hbsag_`cut'==1 & discord_`cut'<.
    local denom = r(N)

    quietly count if hbsag_`cut'==1 & discord_`cut'==0
    local n_conc = r(N)
    local pct_conc = 100*`n_conc'/`denom'

    quietly count if hbsag_`cut'==1 & discord_`cut'==1
    local n_disc = r(N)
    local pct_disc = 100*`n_disc'/`denom'

    post D ("`cut'") ("HBsAg+/anti-HBc+") (`n_conc') (`pct_conc')
    post D ("`cut'") ("HBsAg+/anti-HBc-") (`n_disc') (`pct_disc')
}

postclose D
use `discsum', clear

replace cutoff = "S/CO >1"  if cutoff=="gt1"
replace cutoff = "S/CO >10" if cutoff=="gt10"

format pct %6.1f

capture drop __0*
export excel using "$OUT/discordance/discordance_summary_counts_unweighted.xlsx", ///
    firstrow(variables) replace

capture drop __0*
save "$OUT/discordance/discordance_summary_counts_unweighted.dta", replace

restore

/*---------------------------------------------------------------*
 | 15. HBV DNA among discordant/concordant cases
 |     Descriptive clinical characterisation among cases.
 |     Weighted estimates included where useful, but clinical
 |     medians remain descriptive.
 *---------------------------------------------------------------*/

capture confirm file "$DATA/hbv_cases.dta"

if _rc==0 {

    use "$DERIVED/hbv_cases.dta", clear
    capture drop __0*

    svyset _n [pweight=totwt], singleunit(centered)
    svyset

    capture confirm numeric variable hbv_num
    if _rc {
        capture confirm numeric variable hbsag_result
        if _rc {
            gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9]+(\.[0-9]+)?)")
        }
        else {
            gen double hbv_num = hbsag_result
        }
    }

    cap drop antihbc_clean 
    gen str30 antihbc_clean = lower(trim(anti_hbcii_interpretation))

	cap drop antihbc_pos
    gen byte antihbc_pos = .
    replace antihbc_pos = 1 if antihbc_clean=="reactive"
    replace antihbc_pos = 0 if antihbc_clean=="nonreactive"

    cap drop hbsag_gt1 hbsag_gt10 discord_gt1 discord_gt10

    gen byte hbsag_gt1 = .
    replace hbsag_gt1 = 1 if hbv_num > 1  & hbv_num < .
    replace hbsag_gt1 = 0 if hbv_num <= 1 & hbv_num < .

    gen byte hbsag_gt10 = .
    replace hbsag_gt10 = 1 if hbv_num > 10  & hbv_num < .
    replace hbsag_gt10 = 0 if hbv_num <= 10 & hbv_num < .

    gen byte discord_gt1 = .
    replace discord_gt1 = 1 if hbsag_gt1==1  & antihbc_pos==0
    replace discord_gt1 = 0 if hbsag_gt1==1  & antihbc_pos==1

    gen byte discord_gt10 = .
    replace discord_gt10 = 1 if hbsag_gt10==1 & antihbc_pos==0
    replace discord_gt10 = 0 if hbsag_gt10==1 & antihbc_pos==1

    label define discord_lbl 0 "HBsAg+/anti-HBc+" 1 "HBsAg+/anti-HBc-", replace
    label values discord_gt1 discord_lbl
    label values discord_gt10 discord_lbl

    /* HBV DNA harmonisation */
    capture confirm numeric variable dna_vl
    if _rc {
        gen str20 dna_vl_str = hbv_viral_load_iu_ml
        replace dna_vl_str = "9.9"         if hbv_viral_load_iu_ml=="<10"
        replace dna_vl_str = "1000000001"  if hbv_viral_load_iu_ml==">1000000000"
        replace dna_vl_str = "0"           if hbv_viral_load_iu_ml=="Target Not Detected"
        destring dna_vl_str, gen(dna_vl) force
    }

    cap drop dna_detectable log10_dna sub_dna_gt1 sub_dna_gt10

    gen byte dna_detectable = .
    replace dna_detectable = 0 if dna_vl==0 | dna_vl<10
    replace dna_detectable = 1 if dna_vl>=10 & dna_vl<.

    label define dna_lbl 0 "HBV DNA not detected / <10" 1 "HBV DNA >=10 IU/mL", replace
    label values dna_detectable dna_lbl

    gen double log10_dna = .
    replace log10_dna = log10(dna_vl) if dna_vl>0

    gen byte sub_dna_gt1  = hbsag_gt1==1  & !missing(discord_gt1)  & !missing(dna_detectable)
    gen byte sub_dna_gt10 = hbsag_gt10==1 & !missing(discord_gt10) & !missing(dna_detectable)

    display as text "======================================================"
    display as text "HBV DNA BY DISCORDANCE: S/CO >1"
    display as text "======================================================"
    tab dna_detectable discord_gt1 if hbsag_gt1==1, row col missing
    svy, subpop(sub_dna_gt1): proportion dna_detectable, over(discord_gt1)
    svy, subpop(sub_dna_gt1): logistic dna_detectable i.discord_gt1

    tabstat dna_vl if hbsag_gt1==1, by(discord_gt1) ///
        stat(n mean sd median iqr min max) columns(statistics)

    ranksum dna_vl if hbsag_gt1==1, by(discord_gt1)

    display as text "======================================================"
    display as text "HBV DNA BY DISCORDANCE: S/CO >10"
    display as text "======================================================"
    tab dna_detectable discord_gt10 if hbsag_gt10==1, row col missing
    svy, subpop(sub_dna_gt10): proportion dna_detectable, over(discord_gt10)
    svy, subpop(sub_dna_gt10): logistic dna_detectable i.discord_gt10

    tabstat dna_vl if hbsag_gt10==1, by(discord_gt10) ///
        stat(n mean sd median iqr min max) columns(statistics)

    ranksum dna_vl if hbsag_gt10==1, by(discord_gt10)

    preserve

    keep individualid iintid hbv_num hbsag_gt1 hbsag_gt10 ///
         antihbc_pos discord_gt1 discord_gt10 ///
         hivstat sex hbvdob birthyear onart vl_suppressed ///
         hbv_viral_load_iu_ml dna_vl dna_detectable log10_dna totwt

    keep if hbsag_gt1==1 | hbsag_gt10==1

    capture drop __0*
    export excel using "$OUT/discordance/discordance_hbv_dna_cases.xlsx", ///
        firstrow(variables) replace

    capture drop __0*
    save "$OUT/discordance/discordance_hbv_dna_cases.dta", replace

    restore
}

else {
    display as error "NOTE: $DATA/hbv_cases.dta not found. HBV DNA discordance section was skipped."
}

/*---------------------------------------------------------------*
 | 16. Simple review-table workbook: unweighted observed counts
 |     Weighted estimates are in the Stata log from svy commands.
 *---------------------------------------------------------------*/

use "$DERIVED/vk_evolve.dta", clear
capture drop __0*

svyset _n [pweight=totwt], singleunit(centered)
svyset

capture confirm numeric variable hbv_num
if _rc {
    capture confirm numeric variable hbsag_result
    if _rc {
        gen double hbv_num = real(regexs(1)) if regexm(hbsag_result, "([0-9]+(\.[0-9]+)?)")
    }
    else {
        gen double hbv_num = hbsag_result
    }
}

capture drop antihbc_clean 
gen str30 antihbc_clean = lower(trim(anti_hbcii_interpretation))

gen byte antihbc_pos = .
replace antihbc_pos = 1 if antihbc_clean=="reactive"
replace antihbc_pos = 0 if antihbc_clean=="nonreactive"

gen byte hbsag_gt1  = hbv_num > 1  if hbv_num < .
gen byte hbsag_gt10 = hbv_num > 10 if hbv_num < .

gen byte discord_gt1 = .
replace discord_gt1 = 1 if hbsag_gt1==1  & antihbc_pos==0
replace discord_gt1 = 0 if hbsag_gt1==1  & antihbc_pos==1

gen byte discord_gt10 = .
replace discord_gt10 = 1 if hbsag_gt10==1 & antihbc_pos==0
replace discord_gt10 = 0 if hbsag_gt10==1 & antihbc_pos==1

label values discord_gt1 discord_lbl
label values discord_gt10 discord_lbl

putexcel set "$OUT/discordance/discordance_review_tables.xlsx", replace

putexcel A1 = "Discordance analysis: HBsAg-positive/anti-HBc-negative phenotype"

putexcel A3 = "S/CO >1: unweighted observed counts"
tab discord_gt1 if hbsag_gt1==1, matcell(M1)
putexcel A4 = matrix(M1), names

putexcel A8 = "S/CO >10: unweighted observed counts"
tab discord_gt10 if hbsag_gt10==1, matcell(M10)
putexcel A9 = matrix(M10), names

putexcel A13 = "Weighted estimates are produced in the Stata log using svyset _n [pweight=totwt], singleunit(centered)."

display as result "Discordance analysis complete."
display as result "Outputs saved in: $OUT/discordance"
