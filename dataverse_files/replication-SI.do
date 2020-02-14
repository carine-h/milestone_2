/* Replication for main document for 
	"Does the @realDonaldTrump really matter to financial markets?"
	Allyson Benton and Andrew Philips
	6/26/19
	All models estimated using Stata 15.1
*/
* ------------------------------------------------------------------------------
set scheme burd // you may need to download this via "findit scheme burd"
use "tweetpeso-final.dta", clear
set seed 235095
version 15.1


* ----------- TABLE 1 -------------------------------------------
* Unit root tests suggest % change in peso is stationary
dfuller price_peso, lags(10)
kpss price_peso, not
dfgls price_peso, maxlag(10)
dfgls price_peso, maxlag(10) not

dfuller pctchange_peso, lags(10)
kpss pctchange_peso, not
dfgls pctchange_peso, maxlag(10)
dfgls pctchange_peso, maxlag(10) not
* everything points to I(0)
* ------------------------------------------------------------


* ----------- TABLE 2 -------------------------------------------
* Summary Statistics
gen d_lnusdstock = d.lnusdstock
gen d_MexUS_targetratediff = d.MexUS_targetratediff
* may need to download sutex "findit sutex"
sutex pctchange_peso tweetdum tweet_2cat lnretweet_count lnfavorite_count tweetdum_3to5 tweetdum_8to12 sentimentsyuzhet_neg sentimentsyuzhet_pos pctchangesp500 bondspread10yr_pc d_lnusdstock d_MexUS_targetratediff BdM_any USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother in 3/806
 * ------------------------------------------------------------
 

* -------------FIGURE 1-------------------------------------
* Main independent variable: tweet x log(retweets)
cap drop edate
gen edate = mdy(monthnum, day, year)
format edate %tdnn/dd/YY
tsset edate
twoway line lnretweet_precandidate edate if trump_precandidate == 1, lcolor(black) lwidth(medium) || scatter lnretweet_precandidate edate if trump_precandidate == 1 & lnretweet_precandidate > 0, mcolor(black) msize(medsmall) legend(off) title("Pre-Candidate") xtitle("Time") ytitle("log(Retweet) x Tweet") xlabel(, angle(45))
graph save g1.gph, replace
twoway line lnretweet_candidate edate if trump_candidate == 1, lcolor(black) lwidth(medium) || scatter lnretweet_candidate edate if trump_candidate == 1 & lnretweet_candidate > 0, mcolor(black) msize(medsmall) legend(off) title("Primary Candidate") xtitle("Time") ytitle("log(Retweet) x Tweet") xlabel(, angle(45))
graph save g2.gph, replace
twoway line lnretweet_GOPnominee edate if trump_GOPnominee  == 1, lcolor(black) lwidth(medium) || scatter lnretweet_GOPnominee edate if trump_GOPnominee == 1 & lnretweet_GOPnominee  > 0, mcolor(black) msize(medsmall) legend(off) title("GOP Nominee") xtitle("Time") ytitle("log(Retweet) x Tweet") xlabel(, angle(45))
graph save g3.gph, replace
twoway line lnretweet_lameduck edate if USlameduck  == 1, lcolor(black) lwidth(medium) || scatter lnretweet_lameduck edate if USlameduck == 1 & lnretweet_lameduck > 0, mcolor(black) msize(medsmall) legend(off) title("President-Elect") xtitle("Time") ytitle("log(Retweet) x Tweet") xlabel(, angle(45))
graph save g4.gph, replace
twoway line lnretweet_presidency edate if trump_presidency  == 1, lcolor(black) lwidth(medium) || scatter lnretweet_presidency edate if trump_presidency == 1 & lnretweet_presidency > 0, mcolor(black) msize(medsmall) legend(off) title("President") xtitle("Time") ytitle("log(Retweet) x Tweet") xlabel(, angle(45))
graph save g5.gph, replace
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2) ycommon ysize(5) xsize(10)
graph export "figures/SI-retweet-plot.pdf", as(pdf) replace
* ------------------------------------------------------------

* -------------FIGURE 2-------------------------------------
* Main independent variable: tweet x log(favorites)
twoway line lnfavorite_precandidate edate if trump_precandidate == 1, lcolor(black) lwidth(medium) || scatter lnfavorite_precandidate edate if trump_precandidate == 1 & lnfavorite_precandidate > 0, mcolor(black) msize(medsmall) legend(off) title("Pre-Candidate") xtitle("Time") ytitle("log(Favorite) x Tweet") xlabel(, angle(45))
graph save g1.gph, replace
twoway line lnfavorite_candidate edate if trump_candidate == 1, lcolor(black) lwidth(medium) || scatter lnfavorite_candidate edate if trump_candidate == 1 & lnfavorite_candidate > 0, mcolor(black) msize(medsmall) legend(off) title("Primary Candidate") xtitle("Time") ytitle("log(Favorite) x Tweet") xlabel(, angle(45))
graph save g2.gph, replace
twoway line lnfavorite_GOPnominee edate if trump_GOPnominee  == 1, lcolor(black) lwidth(medium) || scatter lnfavorite_GOPnominee edate if trump_GOPnominee == 1 & lnfavorite_GOPnominee > 0, mcolor(black) msize(medsmall) legend(off) title("GOP Nominee") xtitle("Time") ytitle("log(Favorite) x Tweet") xlabel(, angle(45))
graph save g3.gph, replace
twoway line lnfavorite_lameduck edate if USlameduck  == 1, lcolor(black) lwidth(medium) || scatter lnfavorite_lameduck edate if USlameduck == 1 & lnfavorite_lameduck > 0, mcolor(black) msize(medsmall) legend(off) title("President-Elect") xtitle("Time") ytitle("log(Favorite) x Tweet") xlabel(, angle(45))
graph save g4.gph, replace
twoway line lnfavorite_presidency edate if trump_presidency  == 1, lcolor(black) lwidth(medium) || scatter lnfavorite_presidency edate if trump_presidency == 1 & lnfavorite_presidency > 0, mcolor(black) msize(medsmall) legend(off) title("President") xtitle("Time") ytitle("log(Favorite) x Tweet") xlabel(, angle(45))
graph save g5.gph, replace
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph, rows(2) ycommon ysize(5) xsize(10)
graph export "figures/SI-favorite-plot.pdf", as(pdf) replace
tsset ts_nomiss
* ------------------------------------------------------------



 * ------------------------ FIGURE 3 ------------------------------------
* graphical illustrate of the burn-in period necessary to produce stable pre-shock values
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
* simulate
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
drop time
gen time = _n - 1
* for lame duck
twoway rspike het_lame_ll95_ het_lame_ul95_ time, lcolor(navy) lwidth(medthick) || scatter het_lame_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("") ytitle("Expected Error Variance") legend(off) title("President-Elect") xline(28, lcolor(black) lwidth(thin))
graph save g1.gph, replace
keep in 29/41
drop time
gen time = _n - 1
twoway rspike het_lame_ll95_ het_lame_ul95_ time, lcolor(navy) lwidth(medthick) || scatter het_lame_ time, mcolor(dknavy) msize(large) msymbol(o) xtitle("") ytitle("Expected Error Variance") legend(off) title("President-Elect") xline(28, lcolor(black) lwidth(thin)) 
graph save g2.gph, replace
graph combine g1.gph g2.gph, rows(1) ysize(5) xsize(10)
graph export "figures/SI-burninfigure.pdf", as(pdf) replace
restore
* -----------------------------------------------



* ---------------- FIGURE 4 --------------------------------------------
* Trump's likelihood of winning according to prediction markets (left) and de-trented series (right) 
reg trumppred_WTA ts_nomiss
predict trumppred_WTA_notrend, res
tsset edate
* grab locals for the important dates (primary candidate, GOP nominee)
local primarycand = d(16jun2015) // announces presidential bid (June 16, 2015)
local gopnom = d(19jul2016) // GOP nominee day (JUly 19, 2016)
twoway line trumppred_WTA edate  if trumppred_WTA_notrend != ., ytitle("Pct. Chance of Winning") xtitle("") xlabel(, angle(45)) tlabel(01jan2015 01apr2015 01jul2015  01oct2015    01jan2016 01apr2016 01jul2016  01oct2016  01jan2017) lwidth(medium) xline(`primarycand' `gopnom', lcolor(black) lpattern(solid))
graph save g1.gph, replace
twoway line trumppred_WTA_notrend edate if trumppred_WTA_notrend != ., ytitle("Pct. Chance of Winning (De-Trended)") xtitle("") xlabel(, angle(45)) tlabel(01jan2015 01apr2015 01jul2015  01oct2015 01jan2016 01apr2016 01jul2016  01oct2016  01jan2017) lwidth(medium) xline(`primarycand' `gopnom', lcolor(black) lpattern(solid)) yline(0, lcolor(black) lpattern(solid))
graph save g2.gph, replace
graph combine g1.gph g2.gph, rows(1) ysize(5) xsize(10)
graph export "figures/wta-vs-predictionmkt.pdf", as(pdf) replace
tsset ts_nomiss
* ------------------------------------------------------------

* -------------	TABLE 3 ------------------------------------
* Robustness: prediction markets
dfgls trumppred_WTA
dfgls trumppred_WTA, not  // trend-stationary
arch pctchange_peso l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trumppred_WTA ts_nomiss, arch(1) garch(1) het(tweetdum trumppred_WTA l.pctchangesp500 l.bondspread10yr_pc ts_nomiss) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso tweetdum  tweetdum_WTA l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trumppred_WTA ts_nomiss, arch(1) garch(1) het(tweetdum tweetdum_WTA trumppred_WTA l.pctchangesp500 l.bondspread10yr_pc ts_nomiss) archmlags(1) vce(oim)
* -----------------------------------------------


* -------------	TABLE 4 ------------------------------------
* Additional sentiment analysis
arch pctchange_peso l.pctchange_peso sentimentr_neg_dum_precandidate sentimentr_neg_dum_candidate sentimentr_neg_dum_GOPnominee sentimentr_neg_dum_lameduck sentimentr_neg_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(sentimentr_neg_dum_precandidate sentimentr_neg_dum_candidate sentimentr_neg_dum_GOPnominee sentimentr_neg_dum_lameduck sentimentr_neg_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
 
arch pctchange_peso l.pctchange_peso sentimentr_pos_dum_precandidate sentimentr_pos_dum_candidate sentimentr_pos_dum_GOPnominee sentimentr_pos_dum_lameduck sentimentr_pos_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(sentimentr_pos_dum_precandidate sentimentr_pos_dum_candidate sentimentr_pos_dum_GOPnominee sentimentr_pos_dum_lameduck sentimentr_pos_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim) technique(bfgs)

arch pctchange_peso l.pctchange_peso smtsyuzhet_negcont_precandidate smtsyuzhet_negcont_candidate smtsyuzhet_negcont_GOPnominee smtsyuzhet_negcont_lameduck smtsyuzhet_negcont_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(smtsyuzhet_negcont_precandidate smtsyuzhet_negcont_candidate smtsyuzhet_negcont_GOPnominee smtsyuzhet_negcont_lameduck smtsyuzhet_negcont_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim) technique(bfgs)

arch pctchange_peso l.pctchange_peso smtsyuzhet_poscont_precandidate smtsyuzhet_poscont_candidate smtsyuzhet_poscont_GOPnominee smtsyuzhet_poscont_lameduck smtsyuzhet_poscont_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(smtsyuzhet_poscont_precandidate smtsyuzhet_poscont_candidate smtsyuzhet_poscont_GOPnominee smtsyuzhet_poscont_lameduck smtsyuzhet_poscont_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)

arch pctchange_peso l.pctchange_peso sentimentr_poscont_precandidate sentimentr_poscont_candidate sentimentr_poscont_GOPnominee sentimentr_poscont_lameduck sentimentr_poscont_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(sentimentr_poscont_precandidate sentimentr_poscont_candidate sentimentr_poscont_GOPnominee sentimentr_poscont_lameduck sentimentr_poscont_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) arch0(zero) vce(oim)
* ---------------------------------------------------------------------------


* -------------	TABLE 5 ------------------------------------
* Separating trade- and immigration-specific tweets
arch pctchange_peso l.pctchange_peso traderelev_dum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(traderelev_dum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso traderelev_dum_precandidate traderelev_dum_candidate traderelev_dum_GOPnominee traderelev_dum_lameduck traderelev_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(traderelev_dum_precandidate traderelev_dum_candidate traderelev_dum_GOPnominee traderelev_dum_lameduck traderelev_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero) // Trump didn't tweet about trade during the lame duck period so it SHOULD DROP OUT

arch pctchange_peso l.pctchange_peso immigrantrelev_dum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(immigrantrelev_dum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso immigrantrelev_dum_precandidate immigrantrelev_dum_candidate immigrantrelev_dum_GOPnominee immigrantrelev_dum_lameduck immigrantrelev_dum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(immigrantrelev_dum_precandidate immigrantrelev_dum_candidate immigrantrelev_dum_GOPnominee immigrantrelev_dum_lameduck immigrantrelev_dum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
* -----------------------------------------------



* -------------	TABLE 6 ------------------------------------
* Alternative codings of the peso
arch peso_pctchange_low l.peso_pctchange_low tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
 
arch peso_pctchange_low l.peso_pctchange_low tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero) 

arch peso_pctchange_high l.peso_pctchange_high tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
 
arch peso_pctchange_high l.peso_pctchange_high tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero) 

preserve
sort ts_nomiss
drop if pctchange_FRED == . // this series is a bit shorter than the one we use
gen ts_nomiss2 = _n
tsset ts_nomiss2
arch pctchange_FRED l.pctchange_FRED tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
 
arch pctchange_FRED l.pctchange_FRED tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 USlameduck trump_presidency trump_candidate trump_GOPnominee NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
restore
* -----------------------------------------------


* ----------- TABLE 7 -------------------------------------------
* Number of occurrences across subjects
* To get raw number, multiply frequency by sum. To get number of days in sample, add frequency of != 0's:
* for mexico:
tab mexico
di 34 + 18 + 12 +1 // e.g., 34*1 + 9*2 + 4*3 + 1*4
di 34 + 9 + 4 + 1
* mexico and related:
tab tweet_subject_count
di 87 + (2*67) + (3*22) + (4*28) + (5*9) + (6*5) + (7*7) + (8*3) + 9 + (10*4) + (12*2) + (13*2) + 14 +16
di 87 + 67 + 22 + 28 + 9 + 5 + 7 + 3 + 1 + 4 + 2 + 2 + 1 + 1
* Bernie:
tab tweet_num_bernie
di 25 + (2*7) + (3*2) + 4 + 5 + 6
di 25 + 7 + 2 + 1 + 1 + 1
* Clinton:
tab tweet_num_clinton
di  138  + (2*51) + (3*22) + (4*11) + (5*10) + (6*6) + (7*4) + 8 + (9*3) + (10*5) + (11*2) + 13 + 16 + 28 
di 138 + 51 + 22 + 11 + 10 + 6 + 4 + 1 + 3 + 5 + 2 + 1 + 1 + 1
* Cruz
tab tweet_num_cruz
di 225 + (2*115) + (3*57) + (4*24) + (5*16) + (6*8) + (7*2) + (8*4) + 9 + (10*2)
di 225 + 115 + 57 + 24 + 16 + 8 + 2 + 4 + 1 + 2
* China
tab tweet_num_china
di 36 + (2*9) + (3*2)
di 36 + 9 + 2
* Europe
tab tweet_num_europe
di 9 + (2*1)
di 9 + 1
* ------------------------------------------------------------


* -------------	TABLE 8 ------------------------------------
* placebo subjects
arch pctchange_peso l.pctchange_peso tweet_dum_bernie l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_bernie trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso tweet_dum_clinton l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_clinton trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso tweet_dum_clinton_precandidate tweet_dum_clinton_candidate tweet_dum_clinton_GOPnominee tweet_dum_clinton_lameduck tweet_dum_clinton_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_clinton_precandidate tweet_dum_clinton_candidate tweet_dum_clinton_GOPnominee tweet_dum_clinton_lameduck tweet_dum_clinton_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso tweet_dum_cruz l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_cruz trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch pctchange_peso l.pctchange_peso tweet_dum_cruz_precandidate tweet_dum_cruz_candidate tweet_dum_cruz_GOPnominee tweet_dum_cruz_lameduck tweet_dum_cruz_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_cruz_precandidate tweet_dum_cruz_candidate tweet_dum_cruz_GOPnominee tweet_dum_cruz_lameduck tweet_dum_cruz_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
  
arch pctchange_peso l.pctchange_peso tweet_dum_china l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweet_dum_china trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(xb0wt)
* -----------------------------------------------


* -------------	TABLE 9 ------------------------------------
* Mexican stock market index
reg pctchange_mex_stock l.pctchange_mex_stock l.pctchange_peso tweetdum  l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency  NAFTA_roundsandother

arch pctchange_mex_stock l.pctchange_mex_stock l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) archmlags(1) vce(oim)

arch pctchange_mex_stock l.pctchange_mex_stock l.pctchange_peso tweetdum l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(xb0wt)

arch pctchange_mex_stock l.pctchange_mex_stock l.pctchange_peso tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
* -----------------------------------------------


* -------------	TABLE 10 ------------------------------------
* DV is bond market i-rate spread
reg bondspread10yr_pc l.bondspread10yr_pc l.pctchange_peso tweetdum l.pctchangesp500  d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother

arch bondspread10yr_pc l.bondspread10yr_pc l.pctchange_peso tweetdum  l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) archmlags(1) vce(oim)

arch bondspread10yr_pc l.bondspread10yr_pc l.pctchange_peso tweetdum  l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum trump_candidate trump_GOPnominee USlameduck trump_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)

arch bondspread10yr_pc l.bondspread10yr_pc l.pctchange_peso tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, arch(1) garch(1) het(tweetdum_precandidate tweetdum_candidate tweetdum_GOPnominee tweetdum_lameduck tweetdum_presidency l.pctchangesp500 l.bondspread10yr_pc) archmlags(1) vce(oim) arch0(zero)
* -----------------------------------------------


* -------------	TABLE 11 ------------------------------------
* Asymmetric ARCH/GARCH models
arch pctchange_peso l.pctchange_peso tweetdum l.bondspread10yr_pc l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, saarch(1) arch(1) garch(1) vce(oim) arch0(xb0wt)
generate et = (_n-64)/15
predict sigma2_saarch, variance at(et 0.5)

arch pctchange_peso l.pctchange_peso tweetdum l.bondspread10yr_pc l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, earch(1) egarch(1) vce(oim) arch0(xb0wt)
predict sigma2_egarch, variance at(et 0.5)

arch pctchange_peso l.pctchange_peso tweetdum l.bondspread10yr_pc l.pctchangesp500 d.lnusdstock d.MexUS_targetratediff BdM_any USpres2016 l.USpres2016 trump_candidate trump_GOPnominee USlameduck trump_presidency NAFTA_roundsandother, aparch(1) pgarch(1) vce(oim) arch0(xb0wt)
predict sigma2_aparch, variance at(et 0.5)

* create news response function
twoway line sigma2_saarch et in 3/124, lwidth(vthick) lcolor(black) lpattern(dot) || line sigma2_egarch et in 3/124, lwidth(medthick) lpattern(dash) || line sigma2_aparch et in 3/124, lwidth(medthick) title("") xtitle("Epsilon") ytitle("Expected Conditional Variance") legend(order(1 "SAARCH" 2 "EGARCH" 3 "APARCH"))
graph export "figures/asymmetric-newsresponse.pdf", as(pdf) replace
* -----------------------------------------------


* -------------	TABLE 12 ------------------------------------
* VECM Results show the S&P 500 drives the Mexican stock market
varsoc price_sp500 price_mex_stock
vec price_sp500 price_mex_stock, lags(4) trend(t)
vecstable, graph // looks good
* -----------------------------------------------


* -------------	TABLE 13 ------------------------------------
* Granger causality tests
test [D_price_mex_stock]LD.price_sp500 [D_price_mex_stock]L2D.price_sp500 [D_price_mex_stock]L3D.price_sp500
test [D_price_sp500]LD.price_mex_stock [D_price_sp500]L2D.price_mex_stock [D_price_sp500]L3D.price_mex_stock
* -----------------------------------------------


* -------------	FIGURE 6 ------------------------------------
vec price_sp500 price_mex_stock, lags(4) trend(t)
irf create myirf, set(myirfs) replace step(8)
irf graph irf
irf graph oirf
graph export "figures/vec_oirf.pdf", as(pdf) replace
* -----------------------------------------------


 
