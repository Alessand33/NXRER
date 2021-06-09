* This code is for the reproduction of Figure 6 
* The Dynamics of the U.S. Trade Balance and Real Exchange Rate: The J Curve and Trade Costs?
* By: George Alessandria & Horag Choi

* Need Data_For_Fig6.xlsx for data.
* Step 1: Run Stata program (Figure6_Wedges_Step1.do) to get the data for Step 2
* Step 2: Run GAUSS program (Figure6_Wedges_Step2.pgm) to get the Home and foreign wedges
* Step 3: Run EViews program (Figure6_Wedges_Step3.wf1) for the final commonfactor and its HP filtered data.


clear
cls

import excel "Data_for_Fig6.xlsx", sheet("Sheet1") cellrange(B1:BT237) firstrow

set obs 236
generate time = tq(1957q1) + _n-1
format %tq time
tsset time

gen rern = - rer
gen totqn = - totq

/**********************/
/* ECM3 */
/**********************/
nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totqn - {g_LRNII=-2}*L.NIIY - L.Rel_Inc )+{g_SR=-0.25}*D.totqn +{dniiy2}*d2.NIIY +{dniiy1}*d.NIIY +{IncSR=0.75}*D.Rel_Inc) if time<=tq(2015q4) & time>=tq(1979q1), vce(hac nw)
predict resid_ECM3, residuals,  if time<=tq(2015q4) & time>=tq(1979q1)

gen resid_ECM3srlr = Dnxr - ( -_b[/adj]*(L.nxr -_b[/beta1] -_b[/g_LR]*L.totqn - _b[/g_LRNII]*L.NIIY - L.Rel_Inc) +_b[/g_SR]*D.totqn +_b[/dniiy2]*d2.NIIY +_b[/dniiy1]*d.NIIY +_b[/IncSR]*D.Rel_Inc )
gen resid_ECM3sr  = Dnxr - ( _b[/g_SR]*D.totqn +_b[/dniiy2]*d2.NIIY +_b[/dniiy1]*d.NIIY +_b[/IncSR]*D.Rel_Inc )

gen m 				= ln(M)
gen pc 				= ln(CPIUS)
gen pmpc 			= ln(Pm) - pc
gen logDA 			= ln(DA)

gen resid_home = d.m - ( -_b[/adj]*(L.m -_b[/beta1] +_b[/g_LR]*L.pmpc - _b[/g_LRNII]*L.NIIY - L.logDA) -_b[/g_SR]*D.pmpc +_b[/dniiy2]*d2.NIIY +_b[/dniiy1]*d.NIIY +_b[/IncSR]*D.logDA )
gen resid_foreign = resid_ECM3srlr + resid_home 
