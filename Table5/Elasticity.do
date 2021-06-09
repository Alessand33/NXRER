set more off
clear
cls

import excel "IR_SMPLE_All", sheet("Sheet1") cellrange(A1:E100001) firstrow
*import excel "IR_SMPLE_NoZ", sheet("Sheet1") cellrange(A1:E100001) firstrow
*import excel "IR_SMPLE_NoXi", sheet("Sheet1") cellrange(A1:E100001) firstrow
*import excel "IR_SMPLE_NoB", sheet("Sheet1") cellrange(A1:E100001) firstrow

set obs 100000


generate time = 1 + _n-1
format %ty time
tsset time



	gen EXIMRDSD1 = EXIMR1-DSD1
	gen DEXIMRDSD1 = d.EXIMRDSD1
	gen DEXIMR1 = d.EXIMR1
	gen DTOTQ1 = d.TOTQ1
	gen DDSD1 = d.DSD1
	gen EXIMRYSY1 = EXIMR1-YSY1
	gen DYSY1 = d.YSY1


regress EXIMRDSD1 TOTQ1 if time>=1 & time<=100000
regress d.EXIMRDSD1 d.TOTQ1 if time>=2 & time<=100000
regress DEXIMRDSD1 L.EXIMRDSD1 L.TOTQ1 DTOTQ1  if time>=2 & time<=100000
nl (DEXIMRDSD1 = -{Adj=.1}*(L.EXIMRDSD1-{cons}-{g_LR=5}*L.TOTQ1)+{g_SR=0.25}*DTOTQ1 ) if time>=2 & time<=100000 
