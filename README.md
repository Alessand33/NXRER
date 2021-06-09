# NXRER
Replication notes for 

"The Dynamics of the U.S. Trade Balance and Real Exchange Rate: The J Curve and Trade Costs?" 

By: George Alessandria & Horag Choi

DATA

IPWORLD.dta - IP data for range of US trading partners to extend Dallas Fed IP Measure

DallasFedIP201801.xlsx - IP Data for US and trading partners from Dallas Fed

2016_USDATA.xlsx - US data

BK1.5 - reports cross-correlations between net flows & the real exchange rate for the canonical BKK Model with an Armington Elasticity of 1.5

BK1_5INv - reports cross-correlations between net flows & the real exchange rate for the BKK Model with capital intensive trade with an Armington Elasticity of 1.5

STATA FILES
0_TablesData - produces tables 1 and Appendix Tables 2,3,4 
0_Figure3_BKK - produces figure 3

MODEL FILES

: Benchmark Model 


Checking list:

Root\Models\
   All Models etimated and simulated 

Table 1
   Root\0_TablesData.do 

Table 2
   Root\Table2 
      Table2_Decomposition.txt: Data file for Table 2. Obtained from regressions for Table 1
      Table2_Decomposition.pgm: GAUSS Program file for statistics used reported in Table 2
      Table2_Decomposition.out: output file from Table2_Decomposition.pgm

Table 3
      Root\Models\Benchmark\dynare_Benchmark_results.mat  
      Benchmark parameters & estimates 

Table 4
      Root\Models\
      Alternative model estimates : *_results.mat

Table 5
      Root\Table5 
      Armington Elasticity from the Benchmark Model 
      For details see Root\Table5\Readme_Table5_Steps.txt

Table 6
      Root\Models\* 
      Use the results wih Shock_Decomposition in Conditional variance decompostion
      Take the variance of each shock's contribution and divide them with the data variance.

Figure 1 & 2
      Root\Figure1&2.xlsx
      data from 0_TablesData.do 
      Business cycle data are from the estimation in \Root\Table2\Table2_Decomposition.pgm

Figure 3
      Root\0_Figure3_BKK.do
      Root\Models\BKK\Start_BKK.m to get model moments.	

Figure 4
      Root\Table2\Table2_Decomposition.pgm 
         Panel A.  Output of "NXR and BC components"
         Panel B.  Output of "NXR and TW components"

Figure 5
      Root\Figure5\
         Figure5_TBYDispersionPWT.do

Figure 6
      Root\Figure6\
         Check Readme_Figure6_Steps.txt for detailed steps.

Figure 7
      Root\Figure7\Figure7.m
         Use the results stored in Root\Models\Benchmark\dynare_Benchmark_results.mat

Figure 8
      Root\Figure8
         Use the results stored in Root\Models\Benchmark\dynare_Benchmark_results.mat

Figure 9
      Root\Figure9
         Use the results stored in Root\Models\Benchmark\dynare_Benchmark_results.mat

Figure 10
      Root\Figure10
         Use Benchmark and PTM model results and then they are modified based on the policy choices

Figure 11
      Root\Figure11
         Use the result for Shock Decompostions.

Figure 12
      Root\Figure 12
         Use the results from Start_NoAsymmetric_Estimate.m

Figure 13
      Root\Figure13
         Use the results of Benchmark and NoPTM

Table A1 
      Root\
         Model setup in BKK model

Table A2 
      Root\
         0_TablesData

Table A3 
      Root\0_TablesData.do
Table A4 
       Root\0\TablesData.do

Figure A1 
       Root\Figure A1\**

Figure A2 
       Root\Models\*
         Use the results of varius models

Figure A3
       Root\FigureA3\*.*
         Asymmetric models for the graphs are in this directory


         

estimation(datafile=Data02_world, mode_compute = 6, nobs=139,mh_replic=200000,mh_nblocks=5,mh_jscale=0.3, plot_priors = 1, nodisplay);
