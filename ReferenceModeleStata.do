/* Programme destiné à comparer les resultats avec le programme R */
/* 20/01/2015 : Dédinition du modèle de référence, en relation avec le papier  */


clear
*global root "d:/progs/celine/water"
global root "c:/Chris/progs/celine/water"

cd $root

use "$root\data\FinalFile.dta", replace


global exog "i_under18 log_income i_town i_car b08_locenv_water a2_age i_can i_fra"
global endog "isatis_health"
global instrument "itap_2008 iconcernwatpol_2008"


/* Modèle de référence (Table 4, p 19)  */
sspecialreg i_tap Special, exog($exog) endog($endog) iv($instrument) /*
 */ hetero hetv(i_under18 log_income i_town b08_locenv_water) /*
 */ kdens trim(2.5) 

 
 
/* Modèle de TEST SANS HETERO  */
sspecialreg i_tap Special, exog($exog) endog($endog) iv($instrument) /*
 */ kdens trim(2.5) 

/* Même modèlme mais étape par étapes :  sans hetero  */ 
egen  Mspecial = mean(Special)
replace Special = Special - Mspecial

* etape 1  
regress Special $endog $exog $instrument   /*  <-- Ca colle avec specialre.fitN.R  */ 

predict uhat, resid 
sum uhat   /*  <-- Ca colle avec specialre.fitN.R  */ 

* etape 2
kdens uhat , kernel(epanechnikov) gen(fuhat) at(uhat)

sum fuhat
sum i_tap
* etape 3 : T
gen T = (i_tap - (Special >=0))/fuhat  

sum T

* Trimming de T
loc ptrim = 100 - 2.5
egen tlim = pctile(abs(T)), p(`ptrim')
gen byte ttrm = cond((abs(T) > tlim), ., 1) 

count if ttrm ==.
summ T if ttrm !=.

* etape 4 : 2SLS
ivregress 2sls T $exog ($endog = $instrument) if ttrm !=.

