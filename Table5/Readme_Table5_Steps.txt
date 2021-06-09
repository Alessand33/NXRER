Steps for Table 5 Armington Elasticity Estimation

Run Table5_ArmingtonElasticity.m (Matlab)
    to get artificially generated data for estimations:
         Table5_IR_SMPLE_All.xlsx:   All Shocks
         Table5_IR_SMPLE__NoZ.xlsx:  No Productivity Shocks
         Table5_IR_SMPLE__NoXi.xlsx: No Trade Shocks
         Table5_IR_SMPLE__NoB.xlsx:  No beta Shocks
         
    Need dynare_Benchmark_results.mat obtained from the benchmark Matlab simulation program

Run Elasticity.do (Stata) to get the elasticity estimates
Choose one of the following 4 for each case.
    import excel "Table5_IR_SMPLE_All", sheet("Sheet1") cellrange(A1:E100001) firstrow
    *import excel "Table5_IR_SMPLE_NoZ", sheet("Sheet1") cellrange(A1:E100001) firstrow
    *import excel "Table5_IR_SMPLE_NoXi", sheet("Sheet1") cellrange(A1:E100001) firstrow
    *import excel "Table5_IR_SMPLE_NoB", sheet("Sheet1") cellrange(A1:E100001) firstrow
