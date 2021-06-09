* This code is for the reproduction of Figure A1 
* The Dynamics of the U.S. Trade Balance and Real Exchange Rate: The J Curve and Trade Costs?
* By: George Alessandria & Horag Choi

* Need Data_For_Fig6.xlsx
* Step 1: Run Stata program (FigurA1_Wedges_Step1.do) to get the data for Step 2
* Step 2: Run GAUSS program (FigureA1_Wedges_Step2.pgm) to get the Home and foreign wedges
* Step 3: Run EViews program (FigureA1_Wedges_Step3.wf1) for the final commonfactor.


clear
cls

import excel "Data_for_Fig6.xlsx", sheet("Sheet1") cellrange(B1:BT237) firstrow

set obs 236
generate time = tq(1957q1) + _n-1
format %tq time
tsset time

gen rern = - rer
gen totqn = - totq

gen m 				= ln(M)
gen pc 				= ln(CPIUS)
gen pmpc 			= ln(Pm) - pc
gen logDA 			= ln(DA)



/*
******************** 
   Foreign Demand     
********************
*/

nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.DSD - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.DSD)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY)
predict resid_ECM, residuals, if time<=tq(2015q4) & time>=tq(1979q1)
gen resid_home = d.m - ( -_b[/adj]*(L.m -_b[/beta1] -_b[/g_LR]*L.pmpc - _b[/g_LRNII]*L.NIIY - L.logDA) +_b[/g_SR]*D.pmpc +_b[/dniiy2]*d2.NIIY +_b[/dniiy1]*d.NIIY +_b[/IncSR]*D.logDA )
gen resid_foreing = resid_ECM + resid_home

/*
********************
    No Inventory       
********************
*/
nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc )+{g_SR=-0.25}*D.totq+{IncSR=0.75}*D.Rel_Inc)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY)
predict resid_ECMnoInv, residuals, if time<=tq(2015q4) & time>=tq(1979q1)
gen resid_homenoInv = d.m - ( -_b[/adj]*(L.m -_b[/beta1] -_b[/g_LR]*L.pmpc - L.logDA) +_b[/g_SR]*D.pmpc +_b[/IncSR]*D.logDA )
gen resid_foreingnoInv = resid_ECMnoInv + resid_homenoInv

/*
********************
     Benchmark     
********************
*   /
nl (Dnxr = -{adj=.1}*(L.nxr-{beta1}-{g_LR=-2}*L.totq - L.Rel_Inc - {g_LRNII=-2}*L.NIIY )+{g_SR=-0.25}*D.totq+{dniiy1}*d.NIIY+{dniiy2}*d2.NIIY  +{IncSR=0.75}*D.Rel_Inc)  if time>=tq(1979q1), variables(D.DSD L.DSD D.totq D.Rel_Inc L.nxrworld l.totq D.NIIY D2.NIIY)
predict resid_ECMbench, residuals, if time<=tq(2015q4) & time>=tq(1979q1)
gen resid_homebench = d.m - ( -_b[/adj]*(L.m -_b[/beta1] -_b[/g_LR]*L.pmpc - _b[/g_LRNII]*L.NIIY - L.logDA) +_b[/g_SR]*D.pmpc +_b[/dniiy2]*d2.NIIY +_b[/dniiy1]*d.NIIY +_b[/IncSR]*D.logDA )
gen resid_foreingbench = resid_ECMbench + resid_homebench


