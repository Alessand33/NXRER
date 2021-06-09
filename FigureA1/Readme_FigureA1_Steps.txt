Steps for common Trend Estimations

Run FigureA1_Wedges_Step1.do (Stata) to get home and foreign wedges (in ECM form)
    Need Data_For_Fig6.xlsx

Run FigureA1_Wedges_Step2.pgm (GAUSS) to get home and foreign wedges (in level form)
    Need FigrueA1_Wedges_data.txt obtained from Step 1

Run FigureA1_Wedges_Step3.wf1 (EViews) to get home and common factor
    Need FigureA1_Wedges_Step2.out from Step 2

Final output in EViews for Figure A1 (a):
    wcbenchhf: Common factor for Benchmark
    wcf: Common factor for Foreign Demand
    wcnoinvf:  Common factor for No Inventory

Final output in GAUSS for Figure A1 (b):
    The last printout in FigrueA1_Wedges_Step2.out
    
            
