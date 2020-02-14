/* Replication for main document for 
	"Does the @realDonaldTrump really matter to financial markets?"
	Allyson Benton and Andrew Philips
	6/26/19
	Please note that due to differences in the operating system and types of Stata, results will differ very slightly. This appears to cause different results for those reported in Table 2, Model 2 and Table 3, Model 10, although the log likelihoods remain the same. Users trying to get the exact same output should use Stata SE (version 15.1) for Mac
*/
* ------------------------------------------------------------------------------
set scheme burd // you may need to download this via "findit scheme burd"
use "tweetpeso-final.dta", clear
set seed 235095
version 15.1

* ----------------------
* create residual diagnostic program for ARCH. Enders recommends:
* 1. Standardized Residuals: e_t / h_t (e.g., residual divided by conditional variance). Significant Ljung-Box Q suggests AR in mean equation
* 2. Squared standardized residuals. Significant Ljung-Box Q suggests GARCH effects remaining
cap prog drop archdiag
prog define archdiag
* residual diagnostics: ---
set more on
predict resids, resid // obtain residuals
qui predict cond_var, variance // predict conditional variance
qui gen std_resid = resids/sqrt(cond_var) // enders recommends dividing resid by sqrt(cond_var) to obtain standardized residuals
qui gen std_resid_sq = std_resid^2
di "------------------------------------------------"
di "Standardized Residuals"
*  Ljung-Box Q tests:
wntestq std_resid, lags(1)
wntestq std_resid, lags(2)
wntestq std_resid, lags(3)
more
di ""
di ""
di ""
di "------------------------------------------------"
di "Standardized Residuals, Squared"
* Ljung-Box Q tests:
wntestq std_resid_sq, lags(1)
wntestq std_resid_sq, lags(2)
wntestq std_resid_sq, lags(3)

qui drop resids std_resid cond_var std_resid_sq
set more off
end
* --------


* ---------------------------- FIGURE 1 ---------------------------------------
* plot peso and transformed peso:
gen edate = mdy(monthnum, day, year)
format edate %tdnn/dd/YY
tsset edate
twoway line price_peso edate, title("Peso Exchange Rate", color(black)) ytitle("Pesos per US $", size(medlarge)) xtitle("") xlabel(, angle(45) labsize(medium) ) ylabel(,angle(horizontal)) tlabel(01jan2015 01may2015 01sep2015 01jan2016 01may2016 01sep2016 01jan2017 01may2017 01sep2017 01jan2018, format(%dm-CY)) lwidth(medium) scheme(sjmono) graphregion(color(white))
graph save g1.gph, replace
twoway line pctchange_peso edate, title("% Change Peso", color(black)) ytitle("% Change in Peso per US $", size(medlarge)) xtitle("") xlabel(, angle(45) labsize(medium)) ylabel(,angle(horizontal)) tlabel(01jan2015 01may2015 01sep2015   01jan2016 01may2016 01sep2016 01jan2017 01may2017 01sep2017 01jan2018, format(%dm-CY)) lwidth(medium) lcolor(black) scheme(sjmono) graphregion(color(white))
graph save g2.gph, replace
graph combine g1.gph g2.gph, rows(2)
graph export "figures/peso-dollarexrate.pdf", as(pdf) replace
* ------------------------------------------------------------

* ---------------------------- FIGURE 2a ---------------------------------------
* Subject of Mexico-related tweets
tsset ts_nomiss
graph hbar peso_dum nieto_dum pena_dum deport_dum latino_dum tradedeficit_dum nafta_dum tpp_dum mexican_dum immigrant_dum hispanic_dum trade_dum mexico_dum wall_dum immigration_dum illegal_dum border_dum , ascategory title("Proportion of Days with Tweet Containing:") ytitle("Proportion of Days in Sample", size(medlarge)) yvaroptions( relabel(1 "Peso" 2 "Nieto" 3 "Pena" 4 "Deport" 5 "Latino" 6 "Trade Deficit" 7 "NAFTA" 8 "TPP" 9 "Mexican" 10 "Immigrant" 11 "Hispanic" 12 "Trade" 13 "Mexico" 14 "Wall" 15 "Immigration" 16 "Illegal" 17 "Border") lab(labsize(medium))) ylabel(,labsize(medium))
graph export "figures/barchart_proportions.pdf", as(pdf) replace
* ------------------------------------------------------------

* ---------------------------- FIGURE 2b ---------------------------------------
* Monthly # of tweets w.r.t. Mexico:
preserve
collapse (sum) tweetdum, by(year month monthnum)
sort year monthnum
gen counter = 660
replace counter = counter + (_n-1)
format counter %tm
tsset counter
su tweetdum
local max = `r(max)'
gen shade1 = `max' if counter >= m(2015m6) & counter <= m(2016m7) // primary candidate
gen shade2 = `max' if counter >= m(2016m7) & counter <= m(2016m11) // GOP candidate
gen shade3 = `max' if counter >= m(2016m11) & counter <= m(2017m1) // lame duck
gen shade4 = `max' if counter >= m(2017m1) // as president
twoway area shade1 counter, bcolor(gs7) || area shade2 counter, bcolor(gs13) || area shade3 counter, bcolor(gs4) || area shade4 counter, bcolor(gs9) || line tweetdum counter, ytitle("Monthly Count of Keywords") lwidth(medthick) lcolor(black) xlabel(, angle(45)) tlabel(2015m1 2015m5 2015m9 2016m1 2016m5 2016m9   2017m1 2017m5 2017m9 2018m1) xtitle("") legend(order(1 "Primary Candidate" 2 "GOP Candidate" 3 "President-Elect" 4 "President") rows(1)) 
graph save monthlycount.gph, replace
* NOTE: X-AXIS LABELS ARE MADE BY HAND USING STATA'S GRAPH EDITOR AND THEN SAVED BELOW
*graph export "figures/monthlycount.pdf", as(pdf) replace
restore
* ------------------------------------------------------------



* -------------	TABLE 2 ------------------------------------
* ARCH Effects and GARCH(1,1) Models Using Tweet Dummy
reg pctchange_peso l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother
estat ic
archlm, lags(1/5) 
durbina, lags(1/5) 

arch pctchange_peso l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck  trump_presidency  NAFTA_roundsandother, arch(1) garch(1) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag
* -----------------------------------------------

* ------------------------ FIGURE 3 ------------------------------------
* Predictions from Table 2, Model 4
arch pctchange_peso l.pctchange_peso tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

* grab various means for setx
su l.pctchangesp500, meanonly
loc lpctchangesp500 = r(mean)
su l.bondspread10yr_pc, meanonly
loc lbondspread10yr_pc = r(mean)

preserve
mat B = e(b)
mat VCV = e(V)
local names: colnames VCV
clear
* simulate model variables:
corr2data b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16 b17 b18 b19 b20 b21 b22 b23 b24 b25 b26 b27 b28 b29 b30, n(6000) means(B) cov(VCV)
* Begin simulations
* b21 = HET:tweetdum_precandidate
* b22 = HET: tweetdum_candidate
* b23 = HET: tweetdum_GOPnominee
* b24 = HET: tweetdum_lameduck
* b25 = HET: tweetdum_presidency
* b26 = HET: l.pctchangesp500
* b27 = HET: l.bondspread10yr_pc
* b28 = HET: _cons
* b29 = HET: l.ARCH
* b30 = HET: l.GARCH

* 4 scenarios: tweet during precandidate, candidate, GOPnominee, lameduck and presidency:

* t = 0/30: no tweet, all else at means (assuming epsilon = 0):
gen het_precand_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
gen het_cand_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
gen het_GOPnom_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
gen het_lame_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
gen het_pres_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
forv i = 1/30 {
	foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
	loc pasti = `i' - 1
	su `l'`pasti', meanonly
	gen `l'`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	
	}
}
* at t = 31: tweet:
su het_precand_30, meanonly
gen het_precand_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b21*1) + b30*r(mean) // precandidate
su het_cand_30, meanonly
gen het_cand_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b22*1) + b30*r(mean) // candidate
su het_GOPnom_30, meanonly
gen het_GOPnom_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' +  b23*1) + b30*r(mean) // GOP nom
su het_lame_30, meanonly
gen het_lame_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b24*1) + b30*r(mean) // lame duck
su het_pres_30, meanonly
gen het_pres_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b25*1) + b30*r(mean) // pres
*  when t = 32/50
forv i = 32/50 {
	foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
	loc pasti = `i' - 1
	su `l'`pasti', meanonly
	gen `l'`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	
	}
}
* create UL and LL:
forv i = 0/50 {
	foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
		_pctile `l'`i', p(2.5, 5, 12.5, 87.5, 95, 97.5) 
	gen `l'll95_`i' = r(r1)
	gen `l'll90_`i' = r(r2)
	gen `l'll75_`i' = r(r3)
	gen `l'ul75_`i' = r(r4)
	gen `l'ul90_`i' = r(r5)
	gen `l'ul95_`i' = r(r6)
	}
}
* collapse:
collapse het* //ll* ul*
gen count = _n
reshape long het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ ///
het_precand_ll95_ het_cand_ll95_ het_GOPnom_ll95_ het_lame_ll95_ het_pres_ll95_ ///
het_precand_ll90_ het_cand_ll90_ het_GOPnom_ll90_ het_lame_ll90_ het_pres_ll90_ ///
het_precand_ll75_ het_cand_ll75_ het_GOPnom_ll75_ het_lame_ll75_ het_pres_ll75_ ///
het_precand_ul75_ het_cand_ul75_ het_GOPnom_ul75_ het_lame_ul75_ het_pres_ul75_ ///
het_precand_ul90_ het_cand_ul90_ het_GOPnom_ul90_ het_lame_ul90_ het_pres_ul90_ ///
het_precand_ul95_ het_cand_ul95_ het_GOPnom_ul95_ het_lame_ul95_ het_pres_ul95_, j(time) i(count)
keep in 29/41
drop time
gen time = _n - 1
* For precandidate
twoway rspike het_precand_ll95_ het_precand_ul95_ time, lcolor(navy) lwidth(medthin) || scatter het_precand_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("Day", size(medium)) ytitle("") legend(off) title("Pre-Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(.3(.2).9, labsize(medium))
graph save g1.gph, replace
* for candidate
twoway rspike het_cand_ll95_ het_cand_ul95_ time, lcolor(navy) lwidth(medthin) || scatter het_cand_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("Day", size(medium)) ytitle("") legend(off) title("Primary Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(.3(.2).9, labsize(medium))
graph save g2.gph, replace
* for GOP nominee
twoway rspike het_GOPnom_ll95_ het_GOPnom_ul95_ time, lcolor(navy) lwidth(medthin) || scatter het_GOPnom_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("Day", size(medium)) ytitle("") legend(off) title("GOP Nominee", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(.3(.2).9, labsize(medium))
graph save g3.gph, replace
* for lame duck
twoway rspike het_lame_ll95_ het_lame_ul95_ time, lcolor(navy) lwidth(medthin) || scatter het_lame_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("Day", size(medium)) ytitle("") legend(off) title("President-Elect", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(.3(.2).9, labsize(medium))
graph save g4.gph, replace
* for pres
twoway rspike het_pres_ll95_ het_pres_ul95_ time, lcolor(navy) lwidth(medthin) || scatter het_pres_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("Day", size(medium)) ytitle("") legend(off) title("President", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(.3(.2).9, labsize(medium))
graph save g5.gph, replace
* combine together:
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2) l1(Expected Conditional Error Variance)
graph export "figures/Figure3.pdf", as(pdf) replace
restore
* -----------------------------------------------


* ------------- TABLE 3 ------------------------------------
* Intensity of Tweets, Weighting by Retweets and Favorites
arch pctchange_peso l.pctchange_peso tweet_2cat l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_2cat trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweet_2cat_precandidate tweet_2cat_candidate tweet_2cat_GOPnominee tweet_2cat_lameduck tweet_2cat_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1)  het(tweet_2cat_precandidate tweet_2cat_candidate tweet_2cat_GOPnominee tweet_2cat_lameduck tweet_2cat_presidency l.pctchangesp500  l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag
 
arch pctchange_peso l.pctchange_peso lnretweet_count l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016  trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(lnretweet_count trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso lnretweet_precandidate lnretweet_candidate lnretweet_GOPnominee lnretweet_lameduck lnretweet_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(lnretweet_precandidate lnretweet_candidate lnretweet_GOPnominee lnretweet_lameduck lnretweet_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso lnfavorite_count l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016  trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(lnfavorite_count trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso lnfavorite_precandidate lnfavorite_candidate lnfavorite_GOPnominee lnfavorite_lameduck lnfavorite_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(lnfavorite_precandidate lnfavorite_candidate lnfavorite_GOPnominee lnfavorite_lameduck lnfavorite_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag
* -----------------------------------------------------------


* ----------------- FIGURE 4 -----------------------------------
* Coding the intensity of tweets
arch pctchange_peso l.pctchange_peso tweet_2cat_precandidate tweet_2cat_candidate tweet_2cat_GOPnominee tweet_2cat_lameduck tweet_2cat_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1)  het(tweet_2cat_precandidate tweet_2cat_candidate tweet_2cat_GOPnominee tweet_2cat_lameduck tweet_2cat_presidency l.pctchangesp500  l.bondspread10yr_pc) archmlags(1) vce(oim)  arch0(zero)

* grab various means for setx
su l.pctchangesp500, meanonly
loc lpctchangesp500 = r(mean)
su l.bondspread10yr_pc, meanonly
loc lbondspread10yr_pc = r(mean)

preserve
mat B = e(b)
mat VCV = e(V)
local names: colnames VCV
clear
* simulate variables:
corr2data b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16 b17 b18 b19 b20 b21 b22 b23 b24 b25 b26 b27 b28 b29 b30, n(6000) means(B) cov(VCV)
* Begin simulations
* b21 = HET:tweet_2cat_precandidate
* b22 = HET: tweet_2cat_candidate
* b23 = HET: tweet_2cat_GOPnominee
* b24 = HET: tweet_2cat_lameduck
* b25 = HET: tweet_2cat_presidency
* b26 = HET: l.pctchangesp500
* b27 = HET: l.bondspread10yr_pc
* b28 = HET: _cons
* b29 = HET: l.ARCH
* b30 = HET: l.GARCH
    
* 4 scenarios: tweet during precandidate, candidate, GOPnominee, lameduck and presidency:

* t = 0/30: no tweet, all else at means (assuming epsilon = 0):
foreach z in tw1 tw2 {
	gen het_precand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_cand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_GOPnom_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_lame_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_pres_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
}
forv i = 1/30 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			loc pasti = `i' - 1
			su `l'`z'_`pasti', meanonly
			gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	
		}
	}
}
* at t = 31: tweet:
* For 1 tweet
su het_precand_tw1_30, meanonly
	gen het_precand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*1) + b30*r(mean) // precandidate
	su het_cand_tw1_30, meanonly
	gen het_cand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b22*1) + b30*r(mean) // candidate
	su het_GOPnom_tw1_30, meanonly
	gen het_GOPnom_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b23*1) + b30*r(mean) // GOP nom
	su het_lame_tw1_30, meanonly
	gen het_lame_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b24*1) + b30*r(mean) // lame duck
	su het_pres_tw1_30, meanonly
	gen het_pres_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b25*1) + b30*r(mean) // pres
	
* for 2+ tweets
	su het_precand_tw2_30, meanonly
	gen het_precand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*2) + b30*r(mean) // precandidate
	su het_cand_tw2_30, meanonly
	gen het_cand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b22*2) + b30*r(mean) // candidate
	su het_GOPnom_tw2_30, meanonly
	gen het_GOPnom_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b23*2) + b30*r(mean) // GOP nom
	su het_lame_tw2_30, meanonly
	gen het_lame_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b24*2) + b30*r(mean) // lame duck
	su het_pres_tw2_30, meanonly
	gen het_pres_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b25*2) + b30*r(mean) // pres

*  when t = 32/50
forv i = 32/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
		loc pasti = `i' - 1
		su `l'`z'_`pasti', meanonly
		gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	 
		}
	}
}
* create UL and LL:
forv i = 0/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			_pctile `l'`z'_`i', p(2.5, 5, 12.5, 87.5, 95, 97.5) 
			gen `l'`z'_ll95_`i' = r(r1)
			gen `l'`z'_ll90_`i' = r(r2)
			gen `l'`z'_ll75_`i' = r(r3)
			gen `l'`z'_ul75_`i' = r(r4)
			gen `l'`z'_ul90_`i' = r(r5)
			gen `l'`z'_ul95_`i' = r(r6)
		}
	}
}
* collapse:
collapse het*
gen count = _n 
reshape long het_precand_tw1_ het_cand_tw1_ het_GOPnom_tw1_ het_lame_tw1_ het_pres_tw1_ ///
het_precand_tw1_ll95_ het_cand_tw1_ll95_ het_GOPnom_tw1_ll95_ het_lame_tw1_ll95_ het_pres_tw1_ll95_ ///
het_precand_tw1_ll90_ het_cand_tw1_ll90_ het_GOPnom_tw1_ll90_ het_lame_tw1_ll90_ het_pres_tw1_ll90_ ///
het_precand_tw1_ll75_ het_cand_tw1_ll75_ het_GOPnom_tw1_ll75_ het_lame_tw1_ll75_ het_pres_tw1_ll75_ ///
het_precand_tw1_ul75_ het_cand_tw1_ul75_ het_GOPnom_tw1_ul75_ het_lame_tw1_ul75_ het_pres_tw1_ul75_ ///
het_precand_tw1_ul90_ het_cand_tw1_ul90_ het_GOPnom_tw1_ul90_ het_lame_tw1_ul90_ het_pres_tw1_ul90_ ///
het_precand_tw1_ul95_ het_cand_tw1_ul95_ het_GOPnom_tw1_ul95_ het_lame_tw1_ul95_ het_pres_tw1_ul95_ /// Below here for tw2:
het_precand_tw2_ het_cand_tw2_ het_GOPnom_tw2_ het_lame_tw2_ het_pres_tw2_ ///
het_precand_tw2_ll95_ het_cand_tw2_ll95_ het_GOPnom_tw2_ll95_ het_lame_tw2_ll95_ het_pres_tw2_ll95_ ///
het_precand_tw2_ll90_ het_cand_tw2_ll90_ het_GOPnom_tw2_ll90_ het_lame_tw2_ll90_ het_pres_tw2_ll90_ ///
het_precand_tw2_ll75_ het_cand_tw2_ll75_ het_GOPnom_tw2_ll75_ het_lame_tw2_ll75_ het_pres_tw2_ll75_ ///
het_precand_tw2_ul75_ het_cand_tw2_ul75_ het_GOPnom_tw2_ul75_ het_lame_tw2_ul75_ het_pres_tw2_ul75_ ///
het_precand_tw2_ul90_ het_cand_tw2_ul90_ het_GOPnom_tw2_ul90_ het_lame_tw2_ul90_ het_pres_tw2_ul90_ ///
het_precand_tw2_ul95_ het_cand_tw2_ul95_ het_GOPnom_tw2_ul95_ het_lame_tw2_ul95_ het_pres_tw2_ul95_, j(time) i(count)
keep in 29/41
drop time
gen time = _n - 1
gen time2 = time + .33
* For precandidate
twoway rspike het_precand_tw1_ll95_ het_precand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_precand_tw2_ll95_ het_precand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_precand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_precand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Pre-Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g1.gph, replace
* for candidate
twoway rspike het_cand_tw1_ll95_ het_cand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_cand_tw2_ll95_ het_cand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_cand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_cand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Primary Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g2.gph, replace
* for GOP nominee
twoway rspike het_GOPnom_tw1_ll95_ het_GOPnom_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_GOPnom_tw2_ll95_ het_GOPnom_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_GOPnom_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_GOPnom_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("GOP Nominee", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g3.gph, replace
* for lame duck
twoway rspike het_lame_tw1_ll95_ het_lame_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_lame_tw2_ll95_ het_lame_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_lame_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_lame_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President-Elect", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g4.gph, replace
* for pres
twoway rspike het_pres_tw1_ll95_ het_pres_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_pres_tw2_ll95_ het_pres_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_pres_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_pres_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g5.gph, replace
* combine together:
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2)  l1(Expected Conditional Error Variance) ysize(5) xsize(10) //  b1(Day)
graph export "figures/Figure4.pdf", as(pdf) replace
restore
* -----------------------------------------------


* -------------------FIGURE 5---------------------
* Weighting by ln(retweet)
 arch pctchange_peso l.pctchange_peso lnretweet_precandidate lnretweet_candidate lnretweet_GOPnominee lnretweet_lameduck lnretweet_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(lnretweet_precandidate lnretweet_candidate lnretweet_GOPnominee lnretweet_lameduck lnretweet_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
 
* grab various means for setx
su l.pctchangesp500, meanonly
loc lpctchangesp500 = r(mean)
su l.bondspread10yr_pc, meanonly
loc lbondspread10yr_pc = r(mean)

preserve
mat B = e(b)
mat VCV = e(V)
local names: colnames VCV
clear
* simulate variables:
corr2data b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16 b17 b18 b19 b20 b21 b22 b23 b24 b25 b26 b27 b28 b29 b30, n(6000) means(B) cov(VCV)
* Begin simulations
* b21 = HET:lnretweet_precandidate
* b22 = HET: lnretweet_candidate
* b23 = HET: lnretweet_GOPnominee
* b24 = HET: lnretweet_lameduck
* b25 = HET: lnretweet_presidency
* b26 = HET: l.pctchangesp500
* b27 = HET: l.bondspread10yr_pc
* b28 = HET: _cons
* b29 = HET: l.ARCH
* b30 = HET: l.GARCH
    

* 4 scenarios: tweet during precandidate, candidate, GOPnominee, lameduck and presidency
* Compare the average log retweets (8.2285 ~ log(3746.2102)), compare to 90th percentile: (10.22168 ~ log(27492.816))

* t = 0/30: no tweet, all else at means (assuming epsilon = 0):
foreach z in tw1 tw2 {
	gen het_precand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_cand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_GOPnom_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_lame_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_pres_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
}
forv i = 1/30 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			loc pasti = `i' - 1
			su `l'`z'_`pasti', meanonly
			gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)
		}
	}
}
* at t = 31: tweet:
* For average log retweets
su het_precand_tw1_30, meanonly
	gen het_precand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*8.2285) + b30*r(mean) // precandidate
	su het_cand_tw1_30, meanonly
	gen het_cand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b22*8.2285) + b30*r(mean) // candidate
	su het_GOPnom_tw1_30, meanonly
	gen het_GOPnom_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b23*8.2285) + b30*r(mean) // GOP nom
	su het_lame_tw1_30, meanonly
	gen het_lame_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b24*8.2285) + b30*r(mean) // lame duck
	su het_pres_tw1_30, meanonly
	gen het_pres_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b25*8.2285) + b30*r(mean) // pres
	
* for 90th 'pctile
	su het_precand_tw2_30, meanonly
	gen het_precand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*10.22168) + b30*r(mean) // precandidate
	su het_cand_tw2_30, meanonly
	gen het_cand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b22*10.22168) + b30*r(mean) // candidate
	su het_GOPnom_tw2_30, meanonly
	gen het_GOPnom_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b23*10.22168) + b30*r(mean) // GOP nom
	su het_lame_tw2_30, meanonly
	gen het_lame_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b24*10.22168) + b30*r(mean) // lame duck
	su het_pres_tw2_30, meanonly
	gen het_pres_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b25*10.22168) + b30*r(mean) // pres

*  when t = 32/50
forv i = 32/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
		loc pasti = `i' - 1
		su `l'`z'_`pasti', meanonly
		gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	
		}
	}
}
* create UL and LL:
forv i = 0/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			_pctile `l'`z'_`i', p(2.5, 5, 12.5, 87.5, 95, 97.5) 
			gen `l'`z'_ll95_`i' = r(r1)
			gen `l'`z'_ll90_`i' = r(r2)
			gen `l'`z'_ll75_`i' = r(r3)
			gen `l'`z'_ul75_`i' = r(r4)
			gen `l'`z'_ul90_`i' = r(r5)
			gen `l'`z'_ul95_`i' = r(r6)
		}
	}
}
* collapse:
collapse het*
gen count = _n
reshape long het_precand_tw1_ het_cand_tw1_ het_GOPnom_tw1_ het_lame_tw1_ het_pres_tw1_ ///
het_precand_tw1_ll95_ het_cand_tw1_ll95_ het_GOPnom_tw1_ll95_ het_lame_tw1_ll95_ het_pres_tw1_ll95_ ///
het_precand_tw1_ll90_ het_cand_tw1_ll90_ het_GOPnom_tw1_ll90_ het_lame_tw1_ll90_ het_pres_tw1_ll90_ ///
het_precand_tw1_ll75_ het_cand_tw1_ll75_ het_GOPnom_tw1_ll75_ het_lame_tw1_ll75_ het_pres_tw1_ll75_ ///
het_precand_tw1_ul75_ het_cand_tw1_ul75_ het_GOPnom_tw1_ul75_ het_lame_tw1_ul75_ het_pres_tw1_ul75_ ///
het_precand_tw1_ul90_ het_cand_tw1_ul90_ het_GOPnom_tw1_ul90_ het_lame_tw1_ul90_ het_pres_tw1_ul90_ ///
het_precand_tw1_ul95_ het_cand_tw1_ul95_ het_GOPnom_tw1_ul95_ het_lame_tw1_ul95_ het_pres_tw1_ul95_ /// Below here for tw2:
het_precand_tw2_ het_cand_tw2_ het_GOPnom_tw2_ het_lame_tw2_ het_pres_tw2_ ///
het_precand_tw2_ll95_ het_cand_tw2_ll95_ het_GOPnom_tw2_ll95_ het_lame_tw2_ll95_ het_pres_tw2_ll95_ ///
het_precand_tw2_ll90_ het_cand_tw2_ll90_ het_GOPnom_tw2_ll90_ het_lame_tw2_ll90_ het_pres_tw2_ll90_ ///
het_precand_tw2_ll75_ het_cand_tw2_ll75_ het_GOPnom_tw2_ll75_ het_lame_tw2_ll75_ het_pres_tw2_ll75_ ///
het_precand_tw2_ul75_ het_cand_tw2_ul75_ het_GOPnom_tw2_ul75_ het_lame_tw2_ul75_ het_pres_tw2_ul75_ ///
het_precand_tw2_ul90_ het_cand_tw2_ul90_ het_GOPnom_tw2_ul90_ het_lame_tw2_ul90_ het_pres_tw2_ul90_ ///
het_precand_tw2_ul95_ het_cand_tw2_ul95_ het_GOPnom_tw2_ul95_ het_lame_tw2_ul95_ het_pres_tw2_ul95_, j(time) i(count)
keep in 29/41
drop time
gen time = _n - 1
gen time2 = time + .33
* For precandidate
twoway rspike het_precand_tw1_ll95_ het_precand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_precand_tw2_ll95_ het_precand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_precand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_precand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Pre-Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g1.gph, replace
* for candidate
twoway rspike het_cand_tw1_ll95_ het_cand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_cand_tw2_ll95_ het_cand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_cand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_cand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Primary Candidate", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g2.gph, replace
* for GOP nominee
twoway rspike het_GOPnom_tw1_ll95_ het_GOPnom_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_GOPnom_tw2_ll95_ het_GOPnom_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_GOPnom_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_GOPnom_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("GOP Nominee", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g3.gph, replace
* for lame duck
twoway rspike het_lame_tw1_ll95_ het_lame_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_lame_tw2_ll95_ het_lame_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_lame_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_lame_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President-Elect", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g4.gph, replace
* for pres
twoway rspike het_pres_tw1_ll95_ het_pres_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_pres_tw2_ll95_ het_pres_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_pres_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_pres_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g5.gph, replace
* combine together:
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2)  l1(Expected Conditional Error Variance) ysize(5) xsize(10) //  b1(Day)
graph export "figures/Figure5.pdf", as(pdf) replace
restore
* -----------------------------------------------


* --------------- FIGURE 6 ---------------------------------
* Weighting by ln(favorites)
arch pctchange_peso l.pctchange_peso lnfavorite_precandidate lnfavorite_candidate lnfavorite_GOPnominee lnfavorite_lameduck lnfavorite_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(lnfavorite_precandidate lnfavorite_candidate lnfavorite_GOPnominee lnfavorite_lameduck lnfavorite_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
 
* grab various means for setx
su l.pctchangesp500, meanonly
loc lpctchangesp500 = r(mean)
su l.bondspread10yr_pc, meanonly
loc lbondspread10yr_pc = r(mean)

preserve
mat B = e(b)
mat VCV = e(V)
local names: colnames VCV
clear
* simulate variables:
corr2data b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 b13 b14 b15 b16 b17 b18 b19 b20 b21 b22 b23 b24 b25 b26 b27 b28 b29 b30, n(6000) means(B) cov(VCV)
* Begin simulations
* b21 = HET:lnfavorite_precandidate
* b22 = HET: lnfavorite_candidate
* b23 = HET: lnfavorite_GOPnominee
* b24 = HET: lnfavorite_lameduck
* b25 = HET: lnfavorite_presidency
* b26 = HET: l.pctchangesp500
* b27 = HET: l.bondspread10yr_pc
* b28 = HET: _cons
* b29 = HET: l.ARCH
* b30 = HET: l.GARCH
    

* 4 scenarios: tweet during precandidate, candidate, GOPnominee, lameduck and presidency
* Going to compare the average log favorites (9.196428 ~ log(9861.8396)), compare to 90th percentile: (11.61739 ~ log(111011.6))


* t = 0/30: no tweet, all else at means (assuming epsilon = 0):
foreach z in tw1 tw2 {
	gen het_precand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_cand_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_GOPnom_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_lame_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
	gen het_pres_`z'_0 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500')
}
forv i = 1/30 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			loc pasti = `i' - 1
			su `l'`z'_`pasti', meanonly
			gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)
		}
	}
}
* at t = 31: tweet:
* For average log retweets
su het_precand_tw1_30, meanonly
	gen het_precand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*9.196428) + b30*r(mean) // precandidate
	su het_cand_tw1_30, meanonly
	gen het_cand_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b22*9.196428) + b30*r(mean) // candidate
	su het_GOPnom_tw1_30, meanonly
	gen het_GOPnom_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b23*9.196428) + b30*r(mean) // GOP nom
	su het_lame_tw1_30, meanonly
	gen het_lame_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b24*9.196428) + b30*r(mean) // lame duck
	su het_pres_tw1_30, meanonly
	gen het_pres_tw1_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b25*9.196428) + b30*r(mean) // pres
	
* for 90th 'pctile
	su het_precand_tw2_30, meanonly
	gen het_precand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b21*11.61739) + b30*r(mean) // precandidate
	su het_cand_tw2_30, meanonly
	gen het_cand_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500'  + b22*11.61739) + b30*r(mean) // candidate
	su het_GOPnom_tw2_30, meanonly
	gen het_GOPnom_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b23*11.61739) + b30*r(mean) // GOP nom
	su het_lame_tw2_30, meanonly
	gen het_lame_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b24*11.61739) + b30*r(mean) // lame duck
	su het_pres_tw2_30, meanonly
	gen het_pres_tw2_31 = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500' + b25*11.61739) + b30*r(mean) // pres

*  when t = 32/50
forv i = 32/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
		loc pasti = `i' - 1
		su `l'`z'_`pasti', meanonly
		gen `l'`z'_`i'  = exp(b28 + b27*`lbondspread10yr_pc' + b26*`lpctchangesp500') + b30*r(mean)	
		}
	}
}
* create UL and LL:
forv i = 0/50 {
	foreach z in tw1 tw2 {
		foreach l in het_precand_ het_cand_ het_GOPnom_ het_lame_ het_pres_ {
			_pctile `l'`z'_`i', p(2.5, 5, 12.5, 87.5, 95, 97.5) 
			gen `l'`z'_ll95_`i' = r(r1)
			gen `l'`z'_ll90_`i' = r(r2)
			gen `l'`z'_ll75_`i' = r(r3)
			gen `l'`z'_ul75_`i' = r(r4)
			gen `l'`z'_ul90_`i' = r(r5)
			gen `l'`z'_ul95_`i' = r(r6)
		}
	}
}
* collapse:
collapse het*
gen count = _n
reshape long het_precand_tw1_ het_cand_tw1_ het_GOPnom_tw1_ het_lame_tw1_ het_pres_tw1_ ///
het_precand_tw1_ll95_ het_cand_tw1_ll95_ het_GOPnom_tw1_ll95_ het_lame_tw1_ll95_ het_pres_tw1_ll95_ ///
het_precand_tw1_ll90_ het_cand_tw1_ll90_ het_GOPnom_tw1_ll90_ het_lame_tw1_ll90_ het_pres_tw1_ll90_ ///
het_precand_tw1_ll75_ het_cand_tw1_ll75_ het_GOPnom_tw1_ll75_ het_lame_tw1_ll75_ het_pres_tw1_ll75_ ///
het_precand_tw1_ul75_ het_cand_tw1_ul75_ het_GOPnom_tw1_ul75_ het_lame_tw1_ul75_ het_pres_tw1_ul75_ ///
het_precand_tw1_ul90_ het_cand_tw1_ul90_ het_GOPnom_tw1_ul90_ het_lame_tw1_ul90_ het_pres_tw1_ul90_ ///
het_precand_tw1_ul95_ het_cand_tw1_ul95_ het_GOPnom_tw1_ul95_ het_lame_tw1_ul95_ het_pres_tw1_ul95_ /// Below here for tw2:
het_precand_tw2_ het_cand_tw2_ het_GOPnom_tw2_ het_lame_tw2_ het_pres_tw2_ ///
het_precand_tw2_ll95_ het_cand_tw2_ll95_ het_GOPnom_tw2_ll95_ het_lame_tw2_ll95_ het_pres_tw2_ll95_ ///
het_precand_tw2_ll90_ het_cand_tw2_ll90_ het_GOPnom_tw2_ll90_ het_lame_tw2_ll90_ het_pres_tw2_ll90_ ///
het_precand_tw2_ll75_ het_cand_tw2_ll75_ het_GOPnom_tw2_ll75_ het_lame_tw2_ll75_ het_pres_tw2_ll75_ ///
het_precand_tw2_ul75_ het_cand_tw2_ul75_ het_GOPnom_tw2_ul75_ het_lame_tw2_ul75_ het_pres_tw2_ul75_ ///
het_precand_tw2_ul90_ het_cand_tw2_ul90_ het_GOPnom_tw2_ul90_ het_lame_tw2_ul90_ het_pres_tw2_ul90_ ///
het_precand_tw2_ul95_ het_cand_tw2_ul95_ het_GOPnom_tw2_ul95_ het_lame_tw2_ul95_ het_pres_tw2_ul95_, j(time) i(count)
keep in 29/41
drop time
gen time = _n - 1
gen time2 = time + .33
* For precandidate
twoway rspike het_precand_tw1_ll95_ het_precand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_precand_tw2_ll95_ het_precand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_precand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_precand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Pre-Candidate", position(1) ring(0)) xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g1.gph, replace
* for candidate
twoway rspike het_cand_tw1_ll95_ het_cand_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_cand_tw2_ll95_ het_cand_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_cand_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_cand_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("Primary Candidate", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g2.gph, replace
* for GOP nominee
twoway rspike het_GOPnom_tw1_ll95_ het_GOPnom_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_GOPnom_tw2_ll95_ het_GOPnom_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_GOPnom_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_GOPnom_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("GOP Nominee", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g3.gph, replace
* for lame duck
twoway rspike het_lame_tw1_ll95_ het_lame_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_lame_tw2_ll95_ het_lame_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_lame_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_lame_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President-Elect", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g4.gph, replace
* for pres
twoway rspike het_pres_tw1_ll95_ het_pres_tw1_ul95_ time, lcolor(navy) lwidth(medthin) || rspike het_pres_tw2_ll95_ het_pres_tw2_ul95_ time2, lcolor(dkorange) lwidth(medthin) || scatter het_pres_tw1_ time, mcolor(dknavy) msize(large) msymbol(o) || scatter het_pres_tw2_ time2, mcolor(red) msize(large) msymbol(s) xtitle("Day", size(medium)) ytitle("") legend(off) title("President", position(1) ring(0))  xlabel(0(3)12, labsize(medium)) ylabel(, labsize(medium))
graph save g5.gph, replace
* combine together:
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2)  l1(Expected Conditional Error Variance)  ysize(5) xsize(10) // b1(Day)
graph export "figures/Figure6.pdf", as(pdf) replace

restore
* -----------------------------------------------






* -------------		TABLE 4 -----------------------------------------
* Tone of Tweets
arch pctchange_peso l.pctchange_peso sentimentsyuzhet_neg l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(sentimentsyuzhet_neg trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso smtsyuzhet_neg_dum_precandidate smtsyuzhet_neg_dum_candidate smtsyuzhet_neg_dum_GOPnominee smtsyuzhet_neg_dum_lameduck smtsyuzhet_neg_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(smtsyuzhet_neg_dum_precandidate smtsyuzhet_neg_dum_candidate smtsyuzhet_neg_dum_GOPnominee smtsyuzhet_neg_dum_lameduck smtsyuzhet_neg_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso sentimentsyuzhet_pos l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het( sentimentsyuzhet_pos trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso smtsyuzhet_pos_dum_precandidate smtsyuzhet_pos_dum_candidate smtsyuzhet_pos_dum_GOPnominee smtsyuzhet_pos_dum_lameduck smtsyuzhet_pos_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(smtsyuzhet_pos_dum_precandidate smtsyuzhet_pos_dum_candidate smtsyuzhet_pos_dum_GOPnominee smtsyuzhet_pos_dum_lameduck smtsyuzhet_pos_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
estat ic
archdiag
* ------------------------------------------------------------------

* -------------		TABLE 5 -----------------------------------------
* Active Trading Times
arch pctchange_peso l.pctchange_peso tweetdum_3to5 l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_3to5 trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweetdum_3to5_precandidate tweetdum_3to5_candidate tweetdum_3to5_GOPnominee tweetdum_3to5_lameduck tweetdum_3to5_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_3to5_precandidate tweetdum_3to5_candidate tweetdum_3to5_GOPnominee tweetdum_3to5_lameduck tweetdum_3to5_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweetdum_8to12 l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_8to12 USlameduck trump_presidency trump_candidate trump_GOPnominee l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag

arch pctchange_peso l.pctchange_peso tweetdum_8to12_precandidate tweetdum_8to12_candidate tweetdum_8to12_GOPnominee tweetdum_8to12_lameduck tweetdum_8to12_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_8to12_precandidate tweetdum_8to12_candidate tweetdum_8to12_GOPnominee tweetdum_8to12_lameduck tweetdum_8to12_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
estat ic
archdiag
* -----------------------------------------------------------
