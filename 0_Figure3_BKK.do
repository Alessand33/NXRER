* This code reproduces Figure 3 from
* The Dynamics of the U.S. Trade Balance and Real Exchange Rate: The J Curve and Trade Costs?
* By: George Alessandria & Horag Choi

*cd "D:\Dropbox\Exporting Papers\NX_RER\Data"
*cd "C:\Users\georg\Dropbox\Exporting Papers\NX_RER\Data"
*cd "C:\Users\George Alessandria\Dropbox (Phil Research)\Exporting Papers\NX_RER\Data"
cd "D:\Dropbox (Phil Research)\Exporting Papers\NX_RER\Data"
set more off
clear
cls

/*Data from DallasFed (http://www.dallasfed.org/microsites/institute/dgei/ip.cfm)*/
clear
import excel "DallasFedIP201801.xlsx", sheet("IP") cellrange(G8:K461) firstrow
gen time = qofd(Date)
format time %tq
collapse (mean) WorldexUS AdvancedexUS Emerging US, by(time)
save Dallas2016.dta, replace

clear
import excel "2016_USDATA.xlsx", sheet("Statadata") cellrange(A6:U218) firstrow
foreach x of varlist * {
   label var `x' `"`=`x'[1]'"' 
	}
drop in 1
foreach x of varlist * {
   cap destring `x', replace
	}


gen time = qofd(Date)
format time %tq
replace RERO = 100*RERO

foreach var of varlist RER*{
egen temp = mean(`var') if time==tq(1973q1)
egen temp2 = mean(temp)
replace `var' = 100*`var'/temp2
drop temp*
}
replace RERBOG= RERNBIS if time<=tq(1973q1)

merge time using Dallas2016, sort 
order time
drop Date _merge
/* Now lets extend the world ip data*/
merge time using IPWORLD, sort
rename C110ZI 	IPW
rename C110PC 	CPIW
/*rename C111PC 	CPIUS*/
rename C023ZIS	IPEURO
rename C158ZIS	IPJPN 
rename C156ZIS	IPCAN
rename C112ZIS	IPUK  
rename C273ZI	IPMEX
rename C542ZIS	IPKOREA
rename IPGMFSQ	IPUSMFR
rename IPB50001SQ IPUSTOT	

label var IPUSMFR "Industrial Production: Manufacturing (NAICS), Index 2012=100, Quarterly, Seasonally Adjusted"
label var IPUSTOT "Industrial Production: Total index, Index 2012=100, Quarterly, Seasonally Adjusted"

replace US = IPUSTOT

drop _merge

gen ipww = 0.096*ln(IPMEX) +0.263*ln(IPCAN)+ 0.39*ln(IPJPN) +0.178*ln(IPUK)+0.074*ln(IPKOREA)
gen ipwexus = ln(WorldexUS) 
gen ipadvexus = ln(AdvancedexUS )

foreach var of varlist RERBOG REROECD RERNBIS {
egen temp = mean(`var') if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = ln(`var'/temp2)
drop temp*
}

foreach var of varlist ipw* ipadvexus {
egen temp = mean(`var') if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (`var'-temp2)
drop temp*
}
/*line ipw* time*/

replace ipww = ipwexus if time>=tq(1980q1)
replace ipadvexus  = ipww if time<tq(1980q1)
rename RERBOG rer

gen tot 		= - ln(Pm/Px)
gen totq 		= tot+rer
gen totexoil 	= - ln(Pmexoil/Px)
gen totexoilq 	= totexoil+rer
gen nxr 		= ln(X/M)
gen nxn 		= ln(Px*X/(Pm*M))
gen DA 			= 0.5*ln(CG) + 0.5*ln(ID/25)
gen DAN 		= 3500*(0.5*(CGN/2260) + (IDN/1640))
gen MNDAN 		= (MN-MNPET)/DAN
gen Rel_Inc 	= ipww - DA
gen Rel_IPA 	= ipadvexus - ln(US)
gen Rel_IPW 	= ipww - ln(US)
gen NIIY 		= NII/GDP
gen ipus 		= ln(IPUSMFR)
gen nxrworld 	= nxr - Rel_Inc
gen nxr_ripa 	= nxr - Rel_IPA 
gen nxr_ripw 	= nxr - Rel_IPW 
egen tempp 		= mean(nxrworld)
replace nxrworld =nxrworld - tempp
tsset time
gen Dnxr = d.nxr
gen Dnxrworld = d.nxrworld
estimates clear

foreach var of varlist ipus {
egen temp = mean(`var') if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (`var'-temp2)
drop temp*
}

gen tby = (X-M)/GDP
gen tradey = (X+M)/GDP
gen negtotq = -totq
label var tby "Real Trade Balance (% of GDP)"
label var tradey "Real Trade Share of GDP (X+M)/Y"
graph drop _all

foreach var of varlist tot* Rel* nxr* {
egen temp = mean(`var') if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (`var'-temp2)
drop temp*
}


foreach var of varlist CG IDN {
egen temp = mean(ln(`var')) if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (ln(`var')-temp2)
drop temp*
}


****** LETS EXAMINE THE BIASES FROM USING FOREIGN INDUSTRIAL PRODUCTION (APPENDIX TABLE 4)

gen ip =ln(US)
gen gdp = ln(GDP)
gen yminip = gdp - ip
egen temp = mean((yminip)) if time==tq(1980q1)
egen temp2 = mean(yminip)
replace yminip= ((yminip)-temp2)
drop temp*
gen DSR = ipww + ln(1+tby*exp(yminip)*exp(-Rel_IPW))
gen DSD = DSR - DA
gen DSUS = ip + ln(1-tby*exp(yminip))
gen DSDUS = DSR - DSUS
gen nxrworldDSD = nxr - DSD

*****************************************************************************
*
*			FIGURE 3:
*			REPORT LAG AND LEAD RELATIONSHIP IN THEORY AND DATA
*
*
******************************************************************************


gen rern = - rer
qui xcorr rern nxr if time>=tq(1979q1), table generate(xcrernxr)
qui xcorr rern nxn if time>=tq(1979q1), table generate(xcrernxn)
qui xcorr rern nxr_ripa if time>=tq(1979q1), table  generate(xcrernxr_ripa)
qui xcorr rern nxr_ripw if time>=tq(1979q1), table generate(xcrernxr_ripw)
qui xcorr rern nxrworld if time>=tq(1979q1), table  generate(xcrernxrworlddemand)
qui xcorr rern nxrworldDSD if time>=tq(1979q1), table  generate(xcrernxrworlddemand2)

gen totqn = - totq
qui xcorr totqn nxr if time>=tq(1979q1), table  generate(xctotqnxr)
qui xcorr totqn nxn if time>=tq(1979q1), table  generate(xctotqnxn)
qui xcorr totqn nxr_ripa if time>=tq(1979q1), table  generate(xctotqnxr_ripa)
qui xcorr totqn nxr_ripw if time>=tq(1979q1), table  generate(xctotqnxr_ripw)
qui xcorr totqn nxrworld if time>=tq(1979q1), table  generate(xctotqnxrworlddemand)
qui xcorr totqn nxrworldDSD if time>=tq(1979q1), table generate(xctotqnxrworlddemand2)


keep xc* time
gen lags = -20 
replace lags = l.lags  +1 if time>tq(1957q1) 
drop if lags>20
export excel lags xc* using ".\Xcorrupdated.xls", firstrow(variables) replace 
save Xcorrupdated, replace

graph drop _all
twoway (connected xcrernxr xcrernxrworlddemand lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(red black ) lc(red black ) lpattern( l l - l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) title("Figure 3: Comovement of RER{subscript:t}") name(commovement) legend(order(1 "TB{subscript:t+k}" 2 "TBD*D{subscript:t+k}") ring(0) pos(10) cols(1) region(style(none))) yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) note("Note: Based on US from 1979q1 to 2015q4")
use Xcorrupdated, clear
sort lags
merge lags using BKK1_5, sort
drop _merge
merge lags using BKK1_5Inv, sort
drop _merge

graph drop _all
twoway (connected xcrernxr xmrer xmrerK lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red black ) lc(blue red) lpattern( l - .l l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) title("A. Trade Ratio {subscript:t+k}") name(co_xm) legend(order(1 "Data" 2 "BKK" 3 "BKK Invest. Intensive Trade") ring(0) pos(10) cols(1) region(style(none))) yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) note("Note: Based on US from 1979q1 to 2015q4")
twoway (connected xcrernxrworlddemand /*xmddrer*/ XMDD XMDDINV lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red black ) lc(blue red black ) lpattern( l - .l l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) title("B. Trade-Expenditure Ratio {subscript:t+k}") name(co_xmdd) legend(off) /*legend(order(1 "Data" 2 "Theory") ring(0) pos(10) cols(1) region(style(none)))*/ yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) note("Note: Based on US from 1979q1 to 2015q4")
gr combine co_xm co_xmdd , iscale(*0.8) title("Figure 3: Comovement of RER {subscript:t} with") /*subtitle(US & ROW)*/
gr export "graphics\Fig3_Commovement.wmf", replace
twoway (connected xcrernxr xmrer xmrerK lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red black ) lc(blue red) lpattern( l - .l l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) ytitle(Correlation) /*title("A. Trade Ratio {subscript:t+k}")*/ name(co_xm_notitle) legend(order(1 "Data" 2 "BKK" 3 "BKK Invest. Intensive Trade") ring(0) pos(10) cols(1) region(style(none))) yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) /*note("Note: Based on US from 1979q1 to 2015q4")*/
gr export "graphics\tr_comov.eps", replace
twoway (connected xcrernxrworlddemand /*xmddrer*/ XMDD XMDDINV lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red black ) lc(blue red black ) lpattern( l - .l l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle(Correlation) xtitle(Quarter(k)) /*title("B. Trade-Expenditure Ratio {subscript:t+k}")*/ name(co_xmdd_notitle) legend(off) /*legend(order(1 "Data" 2 "Theory") ring(0) pos(10) cols(1) region(style(none)))*/ yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) /*note("Note: Based on US from 1979q1 to 2015q4")*/
gr export "graphics\trd_comov.eps", replace

* CORRELATION WITY IP to SPENDING ADJUSTMENT
* graph drop _all
* twoway (connected xcrernxr xmrer lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red ) lc(blue red) lpattern( l - l -  )), yscale(r(-0.5 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) title("A. Trade Ratio {subscript:t+k}") name(co_xm) legend(order(1 "Data" 2 "Theory") ring(0) pos(10) cols(1) region(style(none))) yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) note("Note: Based on US from 1979q1 to 2015q4")
* twoway (connected xcrernxrworlddemand xmddrer xcrernxrworlddemand2 lags if lags>=-12 & lags<=12, msymbol(i i i i i ) mcolor(blue red black ) lc(blue red black ) lpattern( l - - l -  )), yscale(r(-0.6 0.9)) xlab(-12(4) 12) ytitle("") xtitle(Quarter(k)) title("B. Trade-Expenditure Ratio {subscript:t+k}") name(co_xmdd) legend(order(1 "Data" 2 "Theory" 3 "Data - adjusted for Foreign Demand") ring(0) pos(10) cols(1) region(style(none))) yline(0, lpattern(dash) lc(black) lw(0.1)) xline(0, lpattern(dash) lc(black) lw(0.1)) note("Note: Based on US from 1979q1 to 2015q4")
* gr combine co_xm co_xmdd , iscale(*0.8) title("Figure 3: Comovement of RER {subscript:t} with") /*subtitle(US & ROW)*/
* gr export "graphics\Fig3_Commovement_withDSD.wmf", replace
