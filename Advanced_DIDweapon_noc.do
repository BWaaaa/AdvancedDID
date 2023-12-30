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

* Conduct a 2-period DID 
use pol, clear
preserve 
	keep if treat_c<=2008
	replace treat_c=1 if treat_c==2008
	collapse (mean) no2 co o3 pm10, by(P treat_c) 
	list, sepby(P)
	xtset P treat_c
	loc text "Pollution by treatment status"
	loc i=1
	foreach v of varlist co no2 pm10 {
	    gen d_`v'=-D.`v'
		tw (con `v' P if treat_c==0) (con `v' P if treat_c==1) (con d_`v' P, lp(longdash)),	sch(s2mono) yt("`v'") xt("Year") xla(2005(1)2018) xli(2008, lp(dash) lc(red*1.5)) yli(0,lp(dash) lc(blue)) leg(lab(1 "Control Zones") lab(2 "Low Emission Zone") lab(3 "difference") r(1)) note("`text'") name("graph`v'", replace) t("`v' across year and treatment")
		graph export "figure`i'_PTH`v'.png", as(png) name("graph`v'") replace
		loc i=`i'+1
	}
	gen d_o3=-D.o3
	tw (con o3 P if treat_c==0) (con o3 P if treat_c==1) (con d_o3 P, lp(longdash)), scheme(s2mono) yt("o3") xt("Year") xla(2005(1)2018) xli(2008, lp(dash) lc(red*1.5)) leg(lab(1 "Control Zones") lab(2 "Low Emission Zone") lab(3 "difference") r(1)) note("`text'") name("grapho3", replace) t("Ozone across year and treatment")
	graph export "figure`i'_PTHo3.png", as(png) name("grapho3") replace
restore  
gen treat=0
replace treat=1 if treat_c==2008
gen Post=0
replace Post=1 if Pe>=2008
gen treat_x_Post=treat*Post
foreach v of varlist co no2 o3 pm10 {
    reg `v' treat_x_Post treat Post if treat_c<=2008, clu(mun)
}

* Conduct a M-period Staggered DID
use allpol, clear
ren date year
xtset station year
enc mun, gen(municipal)
qui reg pm10 stag i.year i.stat, vce(clu mun)
est sto b
esttab b, k(stage1_dummy) b(5) se(5) wide // -1.213
est drop b
xtreg pm10 stag i.year, fe vce(clu mun) // -1.213
tab treat_d, gen(treat_m) 
forvalues k=2(1)11{
    loc i=12-`k'
    ren treat_m`k' F`i'
	lab var F`i' "-`i'" 
}
forvalues k=12(1)22{
    local i=`k'-12
    rename treat_m`k' L`i'
	label var L`i' "`i'" 
}
replace F1=0  
xtreg pm10 F* L* i.year, fe vce(clu mun)
est sto olsdid 
coefplot, keep(F* L*) omit vert coefl(F#=-# L#=#) ciop(recast(rcap)) addplot(li @b @at) xli(10, lp(dash) lc(red*0.4)) yli(0, lp(dash) lc(black*0.4)) yt("pm10") xt("Years to treatment") t("TWFE estimators") sch(s2color) name("twfe", replace)  
drop treat_m* F* L*
graph export "figure5_standardtwfe.png", as(png) name("twfe")

* Bacon decomposition
use pol, clear
xtset station year
ba pm10 stag, vce(clu station) ddetail stub(type) gropt(name(bacd1))
graph export "figure6_bd.png", as(png) name("bacd1")
egen weight = sum(typeS)
replace we  = typeS/we
egen mean_est = sum(typeB*we)
egen group_we = sum(we), by(typecg)
egen mean_gr  = sum(typeB*we/gr), by(typecg)
tw (sc typeB we if typecg==2, mc(gs10)) (li mean_e we if typecg==2, lw(1.0) lc("255 91 54")) (li mean_g we if typecg==2, lw(1.0) lc("00 00 62")), xla(0(0.002)0.004) ysc(r(-3.8,3.8)) yla(-2(2)2) name(evl, replace) t("Earlier vs Later Treated") xt("") yt("Estimates") leg(off)
tw (sc typeB we if typecg==1, mc(gs10)) (li mean_e we if typecg==1, lw(1.0) lc("255 91 54")) (li mean_g we if typecg==1, lw(1.0) lc("00 00 62")), xla(0 0.005) ysc(r(-3.8,3.8)) yla(-2(2)2) name(lve, replace) t("Later vs Earlier Treated") xt("") yt("") leg(off)
tw (sc typeB we if typecg==4, mc(gs10)) (li mean_e we if typecg==4, lw(1.0) lc("255 91 54")) (li mean_g we if typecg==4, lw(1.0) lc("00 00 62")), xla(0.1 0.2) xsc(r(0.04,0.28)) ysc(r(-3.8,3.8)) yla(-2(2)2) name(tvn, replace) t("Treated vs Untreated") xt("") ytitle("") leg( order(2 3) lab(2 "Average Est.") lab(3 "Group Est.") c(1) pos(1) ring(0) s(3))
gr combine evl lve tvn, r(1) note("Weight", place(bottom)) name(combi)
graph export "figure7_bd.png", as(png) name(combi) width(1800) height(800) replace

* CS-robust DID
use Tut5_allpol, clear
ren date year
xtset station year
enc mun, gen(municipal)
set seed 1
csdid pm10, i(stat) t(P) g(treat_c) rs(1) clu(mun) agg(simple) // -1.759
csdid pm10, i(stat) t(P) g(treat_c) rs(1) vce(clu mun) agg(event) 
mat tab_csd = r(table)
estat event, estore(csd)
csdid_plot, yt("pm10") xt("Years to treatment") name("CSDID", replace) sch(s2color) t("CSDID estimators")
graph export "figure8_csdid_pm10.png", as(png) name("CSDID")
mat csd_betas = tab_csd[1..8, 3..23]
gen  zeros  = 0
gen  teners = 11
gen  liner  = _n
egen spaner = seq(), f (-6) to (6) 
coefplot mat(csd_betas), ci((5 6)) vert ciop(recast(rarea) color(gs15)) addplot((li zeros liner, lp(dash) lc(red)) (li spaner teners, lp(dash) lc(black*0.4))) coeflab(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) yla(-6 -3 0 3)  sch(s2color) t("Particulate Matter") name("CSpm10", replace) 
graph export "figure9_csdid_pm10.png", as(png) name("CSpm10")
drop zeros teners liner spaner
csdid co, i(stat) t(P) g(treat_c) rs(1) vce(clu mun) agg(event) 
mat csd_co  = r(table)
csdid no2, i(stat) t(P) g(treat_c) rs(1) vce(clu mun) agg(event) 
mat csd_no2 = r(table) 
csdid o3, i(stat) time(P) g(treat_c) rs(1) vce(clu mun) agg(event) 
mat csd_o3 = r(table) 
gen  zeros   = 0
gen  teners  = 11
gen  liner   = _n
egen spanerc = seq(), f (-1) to (1)
replace spanerc = spanerc/5
egen spanern = seq(), f (-16) to (10)  
egen spanero = seq(), f (-6) to (8) 
egen spanerp = seq(), f (-6) to (6) 
matrix csd_co_betas  = csd_co[1..8, 3..23]
matrix csd_no2_betas = csd_no2[1..8, 3..23]
matrix csd_o3_betas  = csd_o3[1..8, 3..23]
coefplot mat(csd_co_betas), ci((5 6)) vert addplot((li zeros liner, lp(dash) lc(red)) (li spanerc teners, lp(dash) lc(black*0.4))) ciop(recast(rarea) color(gs15)) coeflab(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) yla(-0.2(0.1)0.1) yt("ATT") sch(s2color) t("Carbon Monoxide") name("CSco",replace) 
graph export "figure10_csdid_co.png", as(png) name("CSco")
coefplot mat(csd_no2_betas), ci((5 6)) vert addplot((li zeros liner, lp(dash) lc(red)) (li spanern teners, lp(dash) lc(black*0.4))) ciop(recast(rarea) color(gs15)) coeflab(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) yla(-15(5)10) sch(s2color) t("Nitrogen Dioxide") name("CSno2",replace) 
graph export "figure11_csdid_no2.png", as(png) name("CSno2")
coefplot mat(csd_o3_betas), ci((5 6)) vert addplot((li zeros liner, lp(dash) lc(red)) (li spanero teners, lp(dash) lc(black*0.4))) ciop(recast(rarea) color(gs15)) coeflab(Tm10=-10 Tm9=" " Tm8=" " Tm7=" " Tm6=" " Tm5=-5 Tm4=" " Tm3=" " Tm2=" " Tm1=" " Tp0=0 Tp1=" " Tp2=" " Tp3=" " Tp4=" " Tp5=5 Tp6=" " Tp7=" " Tp8=" " Tp9=" " Tp10=10) yla(-5(2.5)7.5) sch(s2color) t("Ozone") name("CSo3",replace) 
graph export "figure12_csdid_o3.png", as(png) name("CSo3")
gr combine CSco CSno2 CSo3 CSpm10, r(1) note("Years to LEZ implmentation", place(bottom)) name(combi2, replace) 
graph export "figure13_csdid_all.png", as(png) name(combi2) width(1800) height(800) replace

* Other methods DID
** Sun and Abraham method
// ssc install avar
// ssc install reghdfe
// ssc install eventstudyinteract
forvalues k=10(-1)2{ 
    gen g_`k' = treat_dur == -`k'
}
forvalues k=0(1)10{ 
    gen g`k'  = treat_dur == `k'
}
eventstudyinteract pm10 g_* g0-g10, cohort(treat_c) control_cohort(nev) absorb(i.stat i.P) vce(clu mun)
mat beta_sadid=e(b_iw) 
mat vars_sadid=e(V_iw)
mata st_matrix("Ses", sqrt(diagonal(st_matrix("e(V_iw)"))))
mat btab_sadid = beta_sadid\Ses'
mat btab_sadid = btab_sadid[1..2,1..9], [0,0]', btab_sadid[1..2,10..20]
mat coln btab_sadid = "-10" "-9" "-8" "-7" "-6" "-5" "-4" "-3" "-2" "-1" "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"
coefplot mat(btab_sadid[1]), se(btab_sadid[2]) vert ciop(recast(rcap)) sch(s2color) addplot(li @b @at) yli(0, lp(dash) lc(black)) xli(10, lp(dash) lc(red)) xt("Years to treatment") yt("pm10") t("SA-DID method") name("SAdid", replace)
graph export "figure14_SADID.png", as(png) name("SAdid")
** de Chaisemartin and D'Haultfoeuille DID_M estimator: 
did_multiplegt pm10 stat P stag, robust_dynamic dynamic(10) seed(1) cluster(muni) // -1.605
did_multiplegt pm10 stat P stag, robust_dynamic dynamic(10) placebo(6) seed(1) cluster(muni) 
graph export "figure15_ddDIDM.png", as(png) name("Graph")
mat beta_dddid = e(estimates)
mat vars_dddid = e(variances)
** Borusyak et.als method: 
did_imputation pm10 stat P f, cluster(muni) autos delta(1) // -1.84
did_imputation pm10 stat P f, allh pre(7) cluster(muni) autos delta(1) minn(0) 
est sto borusdid
event_plot, ciplottype(rcap) graph_opt(yt("pm10") xt("year from treat") t(Borusyak et al method of DID) xti(-7(1)10) xla(-7(1)10) name(BorusDID, replace) )
graph export "figure16_BorusDID.png", as(png) name("BorusDID")
event_plot borusdid beta_dddid#vars_dddid csd beta_sadid#vars_sadid olsdid, stub_lag(tau# Effect_# Tp# g# L#) stub_lead(pre# Placebo_# Tm# g_# F#) plottype(sc) ciplottype(rcap) together perturb(-0.325(0.13)0.325) trimlead(7) noautolegend graph_opt(t("Event study estimators: pm10", s(medlarge) ) xt("Periods since the event") yt("ATT") xla(-7(1)10) yla(0(1)3) leg(order(1 "Borusyak et al." 3 "de Chaisemartin-D'Haultfoeuille" 5 "Callaway-Sant'Anna" 7 "Sun-Abraham" 9 "TWFE") r(3) region(style(none))) xli(-0.5, lc(gs8) lp(dash)) yli(0, lc(gs8)) name(fivecombine, replace) graphregion(color(white)) bgcolor(white) yla(, angle(horizontal))) lag_opt1(msymbol(O) color(cranberry)) lag_ci_opt1(color(cranberry)) lag_opt2(msymbol(Dh) color(navy)) lag_ci_opt2(color(navy)) lag_opt3(msymbol(Th) color(forest_green)) lag_ci_opt3(color(forest_green)) lag_opt4(msymbol(Sh) color(dkorange)) lag_ci_opt4(color(dkorange)) lag_opt5(msymbol(Oh) color(purple)) lag_ci_opt5(color(purple)) 
graph export "figure17_fivecombine.png", as(png) name("fivecombine")
