Steps for common Trend Estimation

Run Figure6_Wedges_Step1.do (Stata) to get home and foreign wedges (in the ECM form)
    Need Data_For_Fig6.xlsx

Run Figure6_Wedges_Step2.pgm (GAUSS) to get home and foreign wedges (in the level form)
    Need Figrue6_Wedges_data.txt obtained from Step 1

Run Figure6_Wedges_Step3.wf1 ("ecm3factor" in EViews) to get home and common factor
    Need Figure6_Wedges_Step2.out from Step 2

Run HP filter (in EViews) for the estimated common factor.
    Need ECM3C_f in EViews

Final output in EViews for Figure 6:
    ECM3C_f: Estimated Common Factor
    ECM3C_f_trend: HP filtered Common Factor 
    
