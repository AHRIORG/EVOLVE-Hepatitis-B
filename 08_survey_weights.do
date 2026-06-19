************************************************************
* Program Name:  EVOLVE_weights
*
* 1. calculate sampling weights for EVOLVE HepB study
*
* 2. datasets used: 
*
* 3, Datasets created: 
*
* Date created :  29 January 2024
* Date modified:                  
* Author:         K Baisley
* Stata version:  Version 18
* Modified by : L Mazibuko
*************************************************************
**********************************************************************************************

version 18
clear
do "src/stata/_config.do"
tempfile selected_ids sample frame swt
 
capture log close
set mem 100000
set more off
set linesize 240
 
****************************************************************************
****************************************************************************


**==================	
**characteristics of all sampled
**==================	


**ones that are actually sampled 
import excel using "$DATA/list_for_EVOLVE_adult_2004_2005.xlsx", clear firstrow
	
	gen hbvdob=1 if year_cohort=="<1965"
		replace hbvdob=2 if year_cohort=="1965-1974"
         replace hbvdob=3 if year_cohort=="1975-1984"
         replace hbvdob=4 if year_cohort=="1985-1994"
         replace hbvdob=5 if year_cohort=="1995-1999"
         replace hbvdob=6 if year_cohort=="2000-2004"
	label def hbvdob 1"<1965" 2"1965-1974" 3"1975-1984" 4"1985-1994" 5"1995-1999" 6"2000-2004" 7"2005-2009" 8"2010-2014" 9"2015-2020"
 	label val hbvdob hbvdob
	label var hbvdob "Year of birth"
		 
	label def sex 1"Male" 2"Female"
	
	gen sexR=1 if sex=="Male"
		replace sexR=2 if sex=="Female"
	label val sexR sex	
	tab sex sexR
	 
	tab hbvdob sexR
		
	label def hiv_res  1"Positive" 2"Negative"
	gen hivR=1 if hiv_res=="Positive"
		replace hivR=2 if hiv_res=="Negative"
	label val hivR hiv_res
	
	tab hiv_res hivR
**males - 75 from each year & HIV group	
	tab hbvdob hivR if sexR==1
**females - need 75 from each year & HIV group 	
	tab hbvdob hivR if sexR==2

save `selected_ids', replace	

use `selected_ids', clear 
	gen ct=1
	collapse (sum) sample=ct,by(hbvdob hivR sexR)
	label var sample "Number sampled in HIV/age/sex stratum"
	rename sexR sex
	rename hivR hiv_res
save `sample', replace	


**==================
**sampling frame -for design weights (probability of selection)
**==================

**Vukuzazi data
use "$DATA/vukuzazi for sampling.dta", clear
	tab hbvdob hiv_res
**males 
	tab hbvdob hiv_res if sex==1
**females 
	tab hbvdob hiv_res if sex==2
save `frame', replace

use `frame', clear	
	drop if hbvdob==7
	gen ct=1
	collapse (sum) tot=ct,by(hbvdob hiv_res sex)
	label var tot "Total population in HIV/age/sex stratum"
	merge 1:1 hbvdob hiv_res sex using `sample'
	drop _m
	
	gen prob=sample/tot
	
	gen strata=1 if hbvdob==1 & sex==1 & hiv_res==1
		replace strata=2 if hbvdob==1 & sex==1 & hiv_res==2
		replace strata=3 if hbvdob==1 & sex==2 & hiv_res==1
		replace strata=4 if hbvdob==1 & sex==2 & hiv_res==2
local k=5		
forval j=2/6 {		
forval i=1/2 {	
		replace strata=`k' if hbvdob==`j' & sex==`i' & hiv_res==1
local k=`k'+1		
		replace strata=`k' if hbvdob==`j' & sex==`i' & hiv_res==2
local k=`k'+1		
}		
}		
		
	gen dwt=1/prob
	label var dwt "Design weight (prob of selection)"
**checkt the weights - for example, of males born before 1965, a higher % of HIV+ ones were selected than of HIV- (28% vs 7.8%), so downweight the HIV+ (the weight is smaller)
	list hbvdob hiv_res sex sample tot prob dwt, noobs
	
save `swt', replace	


**==================
**non-response weights
**==================

**inverse probability weights from Stephen's work
use "$DATA/Vukuzazi_weights_2022_age_sex.dta", clear
	rename _all, lower
	merge 1:1 individualid using `frame'
	keep if _merge==3
	drop _m
	merge m:1 hbvdob hiv_res sex using `swt'
**these are the younger age group so didn't merge in
	list individual hbvdob hiv_res sex if _m==1, noobs
	drop if _m==1
**I think that we can calculate the total weight as non-response*sampling weight?
      gen totwt=ipw_weights*dwt
	
	tab ipw_weights
**seems to be young men	that have highest weights (so lowest prob of attending)
	tab agecat sex if ipw_weight>3
		tab agecat hiv_res if ipw_weight>3
**highest probability
	tab agecat sex if ipw_weight<1.40

	keep individual ipw_weight dwt totwt
	
save "$DERIVED/EVOLVE_weights.dta", replace		
