**********************************************************
** Title: Stata lab 8-3: Advanced DID methods 
**
** Input:  pol   .dta
**         polall.dta
**
** Output: Figure 1:  Parellel trend test: co
**         Figure 2:  Parellel trend test: no2
**         Figure 3:  Parellel trend test: o3
**         Figure 4:  Parellel trend test: pm10
**         Figure 5:  Standard TWFE OLS event plot: pm10
**         Figure 6:  Bacon decomposition: pm10
**         Figure 7:  Bacon decomposition: pm10 (better)
**         Figure 8:  CSDID estimators: pm10
**         Figure 9:  CSDID estimators: pm10 (better)
**         Figure 10: CSDID estimators: co
**         Figure 11: CSDID estimators: no2
**         Figure 12: CSDID estimators: o3
**         Figure 13: CSDID estimators: combined
**         Figure 14: SADID estimators: pm10
**         Figure 15: DIDM/DIDL estimators: pm10
**         Figure 16: BoruDID estimators: pm10
**         Figure 17: Combined-5 estimators: pm10
**
** Date: 2023/12/15, 12/20, 12/21, 12/22, 12/26, 12/30
** 
** Author: Brian Yisong Wang
**********************************************************

clear 
clear all
set more off 
cd "dir\"

**********************************************************
***       高级计量经济-DID方法(Advanced DID)         *****
***                                                  *****
*** 文档结构：                                       *****
*** L51  - 100 : 基本二乘二DID                       *****
*** L101 - 150 : Staggered DID                       *****
*** L151 - 190 : Bacon Decompsition                  *****
*** L191 - 270 : Callaway Santanna DID               *****
*** L271 - 320 : 其他DID方法(SA, dD, Borusyak)       *****
***              271 Sun Abraham                     *****
***              296 de Chaisemartin D'Haultfoeuille *****
***              309 Borusyak et al                  *****
*** L321 - 340 : 高级DID计量方法整合                 *****
***                                                  *****
**********************************************************

* Conduct a 2-period DID 
use pol, clear

** Draw the common trend graph for the group 2008
preserve 
	keep if treat_cohort <= 2008
	replace treat_cohort = 1 if treat_cohort==2008
	collapse (mean) no2 co o3 pm10 ,  by(Period treat_cohort) 
// 	list, sepby(treat_cohort) 
	list, sepby(Period)
	xtset Period treat_cohort
	local text "Pollution by treatment status"
	local i=1
	foreach var of varlist co no2 pm10 {
	    gen diff_`var' = -D.`var'
		twoway (connected  `var'      Period if treat_cohort==0)  ///
				(connected `var'      Period if treat_cohort==1) ///
				(connected diff_`var' Period, lp(longdash)), ///
				scheme(s2mono) ytitle("`var'") xtitle("Year") ///
				xlabel(2005(1)2018) xline(2008, lp(dash) lc(red*1.5)) ///
				yline(0,lp(dash) lc(blue)) /// 
				legend(label(1 "Control Zones") label(2 "Low Emission Zone") ///
				label(3 "difference") row(1)) note("`text'") name("graph`var'", replace) ///
				title("`var' across year and treatment")
				
		graph export "figure`i'_PTH`var'.png", as(png) name("graph`var'") replace
		local i=`i'+1
	}
	
	gen diff_o3 = -D.o3
	
	twoway (connected o3 Period if treat_cohort==0) (connected o3 Period if treat_cohort==1) (connected diff_o3 Period, lp(longdash)), scheme(s2mono) ytitle("o3") xtitle("Year") xlabel(2005(1)2018) xline(2008, lp(dash) lc(red*1.5)) legend(label(1 "Control Zones") label(2 "Low Emission Zone") label(3 "difference") row(1)) note("`text'") name("grapho3", replace) title("Ozone across year and treatment")
	
	graph export "figure`i'_PTHo3.png", as(png) name("grapho3") replace
	
restore  


** Run the 2-period DID for 2-period ATT
gen treat = 0
replace treat = 1 if treat_cohort==2008

gen Post = 0
replace Post = 1 if Period>=2008

gen treat_x_Post = treat*Post
foreach var of varlist co no2 o3 pm10 {
    reg `var' treat_x_Post treat Post if  treat_cohort<=2008, cluster(mun)
}



* Conduct a M-period Staggered DID
use allpol, clear
rename date year
xtset station year
encode mun, gen(municipal)


** Basic DID regression: for Simple Aggregated M-period ATT
quietly reg pm10 stage1_dummy i.year i.station, vce(cluster mun)
estimates store basicdid
esttab basicdid, keep(stage1_dummy) b(5) se(5) wide // -1.213
// did2s pm10, first_stage(year station) second_stage(stage1_dummy) treatment(stage1_dummy) cluster(station)
estimates drop basicdid


** Using fixed effects is equivalent
xtreg pm10 stage1_dummy i.year, fe vce(cluster mun) // -1.213


** A standard Staggered DID setup: period-by-period TWFE estimates
// Essentially there are more coefficients, more treatment status
// Conduct a standard staggered DID, with -1 as the base period
// yit=b1(y(-3)XT)+b2(y(-2)XT)+b3(y(0)XT)+b4(y(1)XT)+...+m_i+l_t+eit
tab treat_dur, gen(treat_m) // A classicla and tatical approach!
forvalues k=2(1)11{
    local i=12-`k'
    rename treat_m`k' F`i'event
	label var F`i'event "-`i'" // An octotastic way to save some locs!
}
forvalues k=12(1)22{
    local i=`k'-12
    rename treat_m`k' L`i'event
	label var L`i'event "`i'" // An octotastic way to save some locs!
}
replace F1event = 0  // base=-1 is a prevalent benchmark in literature

xtreg pm10 F*event L*event i.year, fe vce(cluster mun)

estimates store olsdid 

coefplot, keep(F*event L*event) omitted vertical coeflabels(F#event=-# L#event=#) ciopts(recast(rcap)) addplot(line @b @at) xline(10, lp(dash) lc(red*0.4)) yline(0, lp(dash) lc(black*0.4)) ytitle("pm10") ylabel(,angle(0) format(%2.1f)) xtitle("Years to treatment") title("TWFE estimators") scheme(s2color) name("twfe", replace)  

drop treat_m* F* L*

graph export "figure5_standardtwfe.png", as(png) name("twfe")



* Show problem: Bacon decomposition: we will show some reversed effects
// same averaged treatment effect across all comparisons
use pol, clear
xtset station year

bacondecomp pm10 stage1_dummy, vce(cluster station) ddetail stub(type) gropt(name(bacd1))

graph export "figure6_bd.png", as(png) name("bacd1")


** to replicate the original graph (as well as a better-looking one)
egen weight        = sum(typeS)
replace weight     = typeS/weight
label var weight   "Bacon decomp Group-wise weight"
gen estimate        = typeB
label var estimate "Bacon decomp Group-wise estimates"

// generate the mean line estimates
egen mean_est       = sum(estimate*weight)
egen group_weight   = sum(weight), by(typecgroup)
egen mean_group     = sum(estimate*weight/group_weight), by(typecgroup)

// compare three groups: early vs late, late vs early, and treated vs never treated. One can see the mojority of weight overwhelmingly stay on treated vs never_treated
twoway (scatter estimate weight if typecgroup==2, mcolor(gs10)) ///
(line mean_est weight if typecgroup==2, lwidth(1.0) lcolor("255 91 54")) ///
(line mean_group weight if typecgroup==2, lwidth(1.0) lcolor("00 00 62")), ///
xlabel(0(0.002)0.004) yscale(range(-3.8,3.8)) ylabel(-2(2)2) name(evl, replace) title("Earlier vs Later Treated") xtitle("") ytitle("Estimates") legend(off)
twoway (scatter estimate weight if typecgroup==1, mcolor(gs10)) ///
(line mean_est weight if typecgroup==1, lwidth(1.0) lcolor("255 91 54")) ///
(line mean_group weight if typecgroup==1, lwidth(1.0) lcolor("00 00 62")), ///
xlabel(0 0.005) yscale(range(-3.8,3.8)) ylabel(-2(2)2) name(lve, replace) title("Later vs Earlier Treated") xtitle("") ytitle("") legend(off)
twoway (scatter estimate weight if typecgroup==4, mcolor(gs10)) ///
(line mean_est weight if typecgroup==4, lwidth(1.0) lcolor("255 91 54")) ///
(line mean_group weight if typecgroup==4, lwidth(1.0) lcolor("00 00 62")), ///
xlabel(0.1 0.2) xscale(range(0.04,0.28)) yscale(range(-3.8,3.8)) ylabel(-2(2)2) name(tvn, replace) title("Treated vs Untreated") xtitle("") ytitle("") legend( order(2 3) label(2 "Average Est.") label(3 "Group Est.") cols(1) position(1) ring(0) size(3))

// combine the above graphs
graph combine evl lve tvn, rows(1) note("Weight", placement(bottom)) name(combi)
graph export "figure7_bd.png", as(png) name(combi) width(1800) height(800) replace
	   
	   

* Staggered DID: CS-robust DID: Love Stata Love R.
// The implementation of CSDID in Stata and R are subtly different in the groups to select when panel is non-perfectly balanced. Therefore, the ATT over all sample is slightly different
// The original paper implements EVERYTHING in R, so we need to revamp in Stata. 
use allpol, clear
rename date year
xtset station year
encode mun, gen(municipal)
set seed 1


** csdid consists of two steps: Step 1 estimate, Step 2 plot. 
// Callaway Sant'Antanna method compares pariwise never-treated with treated_at_T 
// csdid pm10, ivar(station) time(Period) gvar(treat_cohort) rseed(1) agg(simple) cluster(municipal) notyet 
csdid pm10, ivar(station) time(Period) gvar(treat_cohort) rseed(1) cluster(municipal) method(dripw) agg(simple) // Similar ATT avergaed over full sample：-1.759


** For plotting, notice csdid relies on t and t-1 comparisons, so there is an estimator (together with CI) for period -1 and period 0.
csdid pm10, ivar(station) time(Period) gvar(treat_cohort) rseed(1) vce(cluster mun) method(dripw) agg(event) 
matrix tab_csd = r(table) // Fantastic!
// di tab_csd[1,1] // Stata具有一种美！是一种脚本语言的美
estat event, estore(csd) // this produces and stores the estimates at the same time

csdid_plot,  ytitle("pm10") ylabel(,angle(0) format(%2.1f)) xtitle("Years to treatment") name("CSDID", replace) scheme(s2color) title("CSDID estimators")
graph export "figure8_csdid_pm10.png", as(png) name("CSDID")


** To replicate (although estimates differ slightly) the original graph:
matrix csd_betas = tab_csd[1..8, 3..23]
gen  zeros  = 0
gen  teners = 11
gen  liner  = _n
egen spaner = seq(), from (-6) to (6) 

coefplot matrix(csd_betas), ci((5 6)) vertical ciopts(recast(rarea) color(gs15)) addplot((line zeros liner, lp(dash) lc(red)) (line spaner teners, lp(dash) lc(black*0.4))) coeflabels(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) ylabels(-6 -3 0 3)  scheme(s2color) title("Particulate Matter") name("CSpm10",replace) 

graph export "figure9_csdid_pm10.png", as(png) name("CSpm10")

drop zeros teners liner spaner


** Pollutants: CO, NO2, O3 as well 
csdid co, ivar(station) time(Period) gvar(treat_cohort) rseed(1) vce(cluster mun) method(dripw) agg(event) 
matrix csd_co  = r(table)
csdid no2, ivar(station) time(Period) gvar(treat_cohort) rseed(1) vce(cluster mun) method(dripw) agg(event) 
matrix csd_no2 = r(table) 
csdid o3, ivar(station) time(Period) gvar(treat_cohort) rseed(1) vce(cluster mun) method(dripw) agg(event) 
matrix csd_o3 = r(table) 

gen  zeros   = 0
gen  teners  = 11
gen  liner   = _n
egen spanerc = seq(), from (-1) to (1)
replace spanerc = spanerc/5
egen spanern = seq(), from (-16) to (10)  
egen spanero = seq(), from (-6) to (8) 
egen spanerp = seq(), from (-6) to (6) 
matrix csd_co_betas  = csd_co[1..8, 3..23]
matrix csd_no2_betas = csd_no2[1..8, 3..23]
matrix csd_o3_betas  = csd_o3[1..8, 3..23]

coefplot matrix(csd_co_betas), ci((5 6)) vertical addplot((line zeros liner, lp(dash) lc(red)) (line spanerc teners, lp(dash) lc(black*0.4))) ciopts(recast(rarea) color(gs15)) coeflabels(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) ylabels(-0.2(0.1)0.1) ytitle("ATT") scheme(s2color) title("Carbon Monoxide") name("CSco",replace) 

graph export "figure10_csdid_co.png", as(png) name("CSco")

coefplot matrix(csd_no2_betas), ci((5 6)) vertical addplot((line zeros liner, lp(dash) lc(red)) (line spanern teners, lp(dash) lc(black*0.4))) ciopts(recast(rarea) color(gs15)) coeflabels(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) ylabels(-15(5)10) scheme(s2color) title("Nitrogen Dioxide") name("CSno2",replace) 

graph export "figure11_csdid_no2.png", as(png) name("CSno2")

coefplot matrix(csd_o3_betas), ci((5 6)) vertical addplot((line zeros liner, lp(dash) lc(red)) (line spanero teners, lp(dash) lc(black*0.4))) ciopts(recast(rarea) color(gs15)) coeflabels(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) ylabels(-5(2.5)7.5) scheme(s2color) title("Ozone") name("CSo3",replace) 

graph export "figure12_csdid_o3.png", as(png) name("CSo3")

graph combine CSco CSno2 CSo3 CSpm10, rows(1) note("Years to LEZ implmentation", placement(bottom)) name(combi2, replace) 

graph export "figure13_csdid_all.png", as(png) name(combi2) width(1800) height(800) replace



* Other methods DID
** Sun and Abraham method: very similar to CSDID method 
// ssc install avar
// ssc install reghdfe
// ssc install eventstudyinteract
forvalues k=10(-1)2{ 
    gen g_`k' = treat_dur == -`k'
}
forvalues k=0(1)10{ 
    gen g`k'  = treat_dur == `k'
}
eventstudyinteract pm10 g_* g0-g10, cohort(treat_cohort) control_cohort(never_treated) absorb(i.station i.Period) vce(cluster mun)

matrix beta_sadid=e(b_iw) 
matrix vars_sadid=e(V_iw)
mata st_matrix("Ses", sqrt(diagonal(st_matrix("e(V_iw)"))))
matrix btab_sadid = beta_sadid\Ses'
matrix btab_sadid = btab_sadid[1..2,1..9], [0,0]', btab_sadid[1..2,10..20]
matrix colnames btab_sadid = "-10" "-9" "-8" "-7" "-6" "-5" "-4" "-3" "-2" "-1" "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"

coefplot matrix(btab_sadid[1]), se(btab_sadid[2]) vertical ciopts(recast(rcap)) scheme(s2color) addplot(line @b @at) yline(0, lp(dash) lc(black)) xline(10, lp(dash) lc(red)) xtitle("Years to treatment") ytitle("pm10") title("SA-DID method") name("SAdid", replace)

graph export "figure14_SADID.png", as(png) name("SAdid")


** de Chaisemartin and D'Haultfoeuille DID_M estimator: 
// DID_M estimator is comparing the t-1 and t outcome whenever there is a change. And it can be extended to the DID_l estimator where there is a comparison of t-1-l and t period when. Do there is by default DID_0, and DID_1 
// ssc install did_multiplegt
did_multiplegt pm10 station Period stage1_dummy, robust_dynamic dynamic(10) seed(1) cluster(municipal) // Averaged ATT approximate same level -1.605

did_multiplegt pm10 station Period stage1_dummy, robust_dynamic dynamic(10) placebo(6) seed(1) cluster(municipal) 

graph export "figure15_ddDIDM.png", as(png) name("Graph")

matrix beta_dddid = e(estimates)
matrix vars_dddid = e(variances)


** Borusyak et.als method: 
// Borusyak et.al.'s method imputes the potential outcomes of the non-treated and then the treated for Y(0), and then takes an average
did_imputation pm10 station Period first_stage, cluster(municipal) autosample delta(1) // Two-by-two DID is very similar, which is good -1.84

did_imputation pm10 station Period first_stage, allhorizons pretrends(7) cluster(municipal) autosample delta(1) minn(0) // staggered DID. Notice results are sensitive to the number of pretrends. This explains why in any paper they report only partial leads. 

estimates store borusdid

// Plot BorusDID
event_plot, ciplottype(rcap) graph_opt(ytitle("pm10") xtitle("year from treat") title(Borusyak et al method of DID) xtick(-7(1)10) xlabel(-7(1)10) name(BorusDID, replace) )
graph export "figure16_BorusDID.png", as(png) name("BorusDID")



* Merge five estimats (Staggered, CSDID, DIDM, SADID, BoruDID)
est tab _all
event_plot borusdid beta_dddid#vars_dddid csd beta_sadid#vars_sadid olsdid, stub_lag(tau# Effect_# Tp# g# L#event) stub_lead(pre# Placebo_# Tm# g_# F#event) /// 
/// This subbing is intuitive: ensure each # represents the period and are aligned. Others are prefixes and suffixes 
plottype(scatter) ciplottype(rcap) together perturb(-0.325(0.13)0.325) trimlead(7) noautolegend ///
graph_opt(title("Event study estimators: pm10", size(medlarge) ) ///
xtitle("Periods since the event") ytitle("ATT") xlabel(-7(1)10) ylabel(0(1)3) ///
legend(order(1 "Borusyak et al." 3 "de Chaisemartin-D'Haultfoeuille" ///
	5 "Callaway-Sant'Anna" 7 "Sun-Abraham" 9 "TWFE") rows(3) region(style(none))) ///
xline(-0.5, lc(gs8) lp(dash)) yline(0, lcolor(gs8)) name(fivecombine, replace) /// 
graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal))) ///
	lag_opt1(msymbol(O) color(cranberry)) lag_ci_opt1(color(cranberry)) ///
	lag_opt2(msymbol(Dh) color(navy)) lag_ci_opt2(color(navy)) ///
	lag_opt3(msymbol(Th) color(forest_green)) lag_ci_opt3(color(forest_green)) ///
	lag_opt4(msymbol(Sh) color(dkorange)) lag_ci_opt4(color(dkorange)) ///
	lag_opt5(msymbol(Oh) color(purple)) lag_ci_opt5(color(purple)) 

graph export "figure17_fivecombine.png", as(png) name("fivecombine")

