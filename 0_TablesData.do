* This code reproduces the Data Tables from
* The Dynamics of the U.S. Trade Balance and Real Exchange Rate: The J Curve and Trade Costs?
* By: George Alessandria & Horag Choi

* Table 1 - Estimates of US Export-Import Ratio
* Table A2 - Estimates of US Export-Import Ratio with Alternative Measures of US Demand (SR)
* Tabl2 A3 - Estimates of US Export-Import Ratio with Alternative Measures of US Demand (SR/LR)
* Tabl2 A4 - Estimates of US Export-Import Ratio adjusting for Spending

*cd "D:\Dropbox\Exporting Papers\NX_RER\Data"
*cd "C:\Users\georg\Dropbox\Exporting Papers\NX_RER\Data"
*cd "C:\Users\George Alessandria\Dropbox (Phil Research)\Exporting Papers\NX_RER\Data"
cd "D:\Dropbox (Phil Research)\Exporting Papers\NX_RER\Replication"

set more off
clear
cls

*use save AC_NXRER_DATA, clear

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

/*
* IMPORT SOME PRICES TO ADD A TERM WITH THE RELATIVE PRICE OF THE WEIGHTED FINAL GOOD to CONSUMPTION

import fred A006RD3Q086SBEA A007RD3Q086SBEA DGDSRD3Q086SBEA DPCERD3Q086SBEA PCDG, clear
gen PGPC = ln(DGDSRD3Q086SBEA/ DPCERD3Q086SBEA)
gen PIPC = ln(A006RD3Q086SBEA / DPCERD3Q086SBEA)
gen PIFPC = ln(A007RD3Q086SBEA / DPCERD3Q086SBEA)
gen time = qofd(daten)
format time %tq PCDG
keep P* time
sort time
save RelativePrices.dta, replace*/

sort time
merge time using RelativePrices.dta
sort time
tsset time
rename PCDG CDG

gen tby = (X-M)/GDP
gen tradey = (X+M)/GDP
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
gen DADUR 		= 0.5*ln(CDG) + 0.5*ln(ID/25)
gen PD 			= 0.5*(PGPC+PIPC)

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
/*line totq nxrworld time*/
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

foreach var of varlist tot* Rel* nxr* DA{
egen temp = mean(`var') if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (`var'-temp2)
drop temp*
}


foreach var of varlist CG IDN ID Ifixed CDG{
egen temp = mean(ln(`var')) if time==tq(1980q1)
egen temp2 = mean(temp)
replace `var' = (ln(`var')-temp2)
drop temp*
}
gen Rel_IncCG = ipww - CG
gen Rel_IncID = ipww - ID
gen Rel_IncIF = ipww - 0.5*(ID+CG)
gen Rel_IncIDUR = ipww - DADUR

*replace DA = ID + CG
gen Rel_IncU = ipww - DA
*drop if time<tq(1979q1) 
*drop if time>tq(2015q4)



*****************************************

* TABLE 1

*****************************************

estimates clear
qui nl (nxr = {beta1}-3*totq + Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4Lmin3
qui nl (nxr = {beta1}-3*tot+ Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4Lmin3tot
qui nl (nxr = {beta1}-3*rer + Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4Lmin3rer
qui nl (nxr = {beta1}+{g_SR=-0.25}*totq + Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4L
qui nl (nxr = {beta1}+{g_SR=-0.25}*tot + Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4Lt
qui nl (nxr = {beta1}+{g_SR=-0.25}*rer + Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4Lq
qui nl (nxr = {beta1}-3*totq + {Inc_LR}*Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(totq Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store ns2015q4LINC
qui nl (Dnxr = {beta1}+{g_SR=-0.25}*D.totq +D.Rel_Inc ) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store nn1post79nws
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq +D.Rel_Inc ) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store nn1post79nw
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq +D.Rel_Inc +{dniiy1}*d.NIIY) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store nn1post79nwni1 
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq +D.Rel_Inc +{dniiy1=1}*d.NIIY+{dniiy2=1}*d2.NIIY) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq d2.NIIY) vce(hac nw)
estimates store nn1post79nwni2
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq +{IncSR=0.75}*D.Rel_Inc) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq) vce(hac nw)
estimates store nn1post79n
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc) if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq d.NIIY) vce(hac nw)
estimates store nn1post79
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {IncLR=1}*L.Rel_Inc)+{g_SR=-0.25}*D.totq +{dniiy1}*d.NIIY+{IncSR=0.75}*D.Rel_Inc)  if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY) vce(hac nw)
estimates store nn2post79
qui nl (Dnxr = -{adj=.1}*(L.nxrworld-{beta1}-{g_LR=-2}*L.totq)+{g_SR=-0.25}*D.totq +{dniiy2=1}*d2.NIIY +{dniiy1=1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY) vce(hac nw)
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {IncLR=1}*L.Rel_Inc - {g_LRNII=-2}*L.NIIY)+{g_SR=-0.25}*D.totq+{dniiy2}*d2.NIIY  +{dniiy1}*d.NIIY+{IncSR=0.75}*D.Rel_Inc)  if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc - {g_LRNII=-2}*L.NIIY)+{g_SR=-0.25}*D.totq+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.Rel_Inc)  if time>=tq(1979q1), variables(D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l2

************************** THIS IS THE BEST REGRESSION
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {g_LRNII=-2}*L.NIIY - L.Rel_Inc         )+{g_SR=-0.25}*D.totq +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY)  vce(hac nw)
estimates store nn1_post79LNI
************************** THIS IS THE BEST REGRESSION
esttab ns2015q4Lmin3 ns2015q4L  nn1post79nws nn1post79nw nn1post79n nn1_post79LNI using "Tables\TabECM_finalNW.tex", nonum ti(Error Correction Model with extended window) se star stats(N rmse r2_a) noeqli mtitles(Level Differences ECM1 ECM2 ECM3)  compress replace
esttab ns2015q4Lmin3 ns2015q4L  nn1post79nws nn1post79nw nn1post79n nn1_post79LNI using "Tables\TabECM_finalNW.csv", nonum ti(Error Correction Model with extended window) se star stats(N rmse r2_a) noeqli mtitles(Level Differences ECM1 ECM2 ECM3)  compress replace
esttab ns2015q4Lmin3 ns2015q4L  nn1post79nws nn1post79nw nn1post79n nn1_post79LNI , t(%7.3f) star stats(N r2_a rmse) noeqli nonum compress
esttab nn1_post79LNI nn2post79l , se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress

**********************************

*			APPENDIX TABLE 4 - EXAMINING THE BIASES FROM USING FOREIGN INDUSTRIAL PRODUCTION

********************************************
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

estimates clear
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {g_LRNII=-2}*L.NIIY - L.Rel_Inc         )+{g_SR=-0.25}*D.totq +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY) vce(hac nw)
estimates store nn1_post79LNI
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {g_LRNII=-2}*L.NIIY - {INCLR}*L.Rel_Inc         )+{g_SR=-0.25}*D.totq +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY) vce(hac nw)
estimates store nn1_post79LNILR
*qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {g_LRNII=-2}*L.NIIY - {INCLR}*L.Rel_IPW )+{g_SR=-0.25}*D.totq +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_IPW)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY)
*estimates store nn1_post79LNIRelIP
*qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {g_LRNII=-2}*L.NIIY - L.Rel_IPW )+{g_SR=-0.25}*D.totq +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_IPW)if time>=tq(1979q1), variables(D.totq l.totq l.nxrworld d.NIIY d2.NIIY)
*estimates store nn1_post79LNIRelIPun
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.DSD - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSD)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l2DSD
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {INCLR}*L.DSD - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSD)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l2DSDLR
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.DSD  )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSD)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l2DSDnoinvLR
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {INCLR}*L.DSD  )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSD)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY) vce(hac nw)
estimates store nn2post79l2DSDLRnoinvLR
*nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.DSDUS - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSDUS)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY)
*estimates store nn2post79l2DSDUS
*nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - {INCLR}*L.DSDUS - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSDUS)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY)
*estimates store nn2post79l2DSDUSLR
esttab nn* , t(%7.2f) star stats(N r2_a rmse) noeqli nonum compress
esttab nn* using "Tables\TabECMDATAAPPENDIX_finalNW.csv", nonum ti(Error Correction Model with extended window) se star stats(N rmse r2_a) noeqli mtitles(IPROW IPROW DEMAND1 DEMAND2 DEMAND3 DEMAND4)  compress replace


********************************

* Table A2:

**************************************

estimates drop _all
qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*D.Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc) vce(hac nw)
estimates store diff
qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*D.Rel_IncCG) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc) vce(hac nw)
estimates store diffCG
qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*D.Rel_IncID) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncID) vce(hac nw)
estimates store diffID
qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*D.Rel_IncIDUR) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncIDUR) vce(hac nw)
estimates store diffIDUR
*qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*D.Rel_IncIF) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncIF)
*estimates store diffIF
qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*(d.ipww - ({alp = 0.5}*D.ID + (1-{alp = 0.5})*D.CG))) if time<=tq(2015q4)& time>=tq(1979q1) , variables(D.totq D.ID) vce(hac nw)
estimates store diffIalp
*qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*(d.ipww - ({alp = 0.5}*D.If + (1-{alp = 0.5})*D.CG))) if time<=tq(2015q4)& time>=tq(1979q1) , variables(D.totq D.If)
*estimates store diffIwNII
*qui nl (Dnxr = {g_SR=-0.25}*D.totq +{gy_SR=-0.25}*(d.ipww - D.ID) ) if time<=tq(2015q4)& time>=tq(1979q1) , variables(D.totq D.ID)
*estimates store diffIw
esttab diff* , se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress mtitles(Benchmark C_Goods Investment CD+ID CgID)
esttab diff* using "Tables\TabDiff_Appendix2.tex", nonum ti(Short-run Income and Price Estimates with Alternative Measures of US Final Demand) se star stats(N rmse r2_a) noeqli mtitles(Benchmark C_Goods Investment CD+ID CgID)  compress replace
esttab diff* using "Tables\TabDiff_Appendix2.csv", nonum ti(Short-run Income and Price Estimates with Alternative Measures of US Final Demand) se star stats(N rmse r2_a) noeqli mtitles(Benchmark C_Goods Investment CD+ID CgID)  compress replace

********************************

* Table A3:

**************************************
estimates drop _all
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc) + {g_SR=-0.25}*D.totq +D.Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc) vce(hac nw)
estimates store ECM
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc) + {g_SR=-0.25}*D.totq +{IncSR=0.75}*D.Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc) vce(hac nw)
estimates store ECMINCSR
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*D.Rel_Inc -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc d2.NIIY) vce(hac nw)
estimates store ECMNIIY
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_Inc - {g_LRNII=-2}*L.NIIY) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_Inc -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc d2.NIIY) vce(hac nw)
estimates store ECM2NIIYLR
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_IncCG - {g_LRNII=-2}*L.NIIY) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_IncCG -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncCG d2.NIIY) vce(hac nw)
estimates store ECM2NIIYLR_CG
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_IncID - {g_LRNII=-2}*L.NIIY) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_IncID -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncID d2.NIIY) vce(hac nw)
estimates store ECM2NIIYLR_ID
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_IncIDUR - {g_LRNII=-2}*L.NIIY) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_IncIDUR -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncCG d2.NIIY) vce(hac nw)
estimates store ECM2NIIYLR_DUR
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_IncCG ) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_IncCG ) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncCG d2.NIIY) vce(hac nw)
estimates store ECM2_CG
qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_IncIDUR ) + {g_SR=-0.25}*D.totq + {IncSR=0.75}*d.Rel_IncIDUR ) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_IncCG d2.NIIY) vce(hac nw)
estimates store ECM2_DUR
esttab ECM* , se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress mtitles(ECM1 ECM1INCSR ECM1INCSR2 ECM1INCSR2NIILR Investment Cons Cons_exNII IPW)
esttab ECM* , se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress
esttab ECM2NIIYLR* , se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress

qui nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-4}*L.totq - L.Rel_Inc - {g_LRNII=-2}*L.NIIY) + {g_SR=-0.25}*D.totq + {g_pgSR1=-0.25}*D.PD +{g_pgSR2=-0.25}*D2.PD +{IncSR=0.75}*d.Rel_Inc -{dniiy1}*d.NIIY -{dniiy2}*d2.NIIY) if time<=tq(2015q4) & time>=tq(1979q1), variables(D.totq D.Rel_Inc d2.NIIY d.PD d2.PD) vce(hac nw)
estimates store ECM2Relp2

esttab ECM2NIIYLR* ECM2Relp2, se(%7.3f) star stats(N r2_a rmse) noeqli nonum compress
esttab ECM2NIIYLR* ECM2Relp2 using "Tables\TabECM_Appendix2.tex", nonum ti(Error Correction Estimates with Alternative Measures of US Final Demand) se star stats(N rmse r2_a) noeqli mtitles(Benchmark C_Goods Investment CD+ID RelPrices)  compress replace
esttab ECM2NIIYLR* ECM2Relp2 using "Tables\TabECM_Appendix2.csv", nonum ti(Error Correction Estimates with Alternative Measures of US Final Demand) se star stats(N rmse r2_a) noeqli mtitles(Benchmark C_Goods Investment CD+ID RelPrices)  compress replace
