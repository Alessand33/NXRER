
clear all;
close all;
clc;

x=load('dynare_NoPTM_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)

Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 
Decomp_DSD=permute(x.oo_.shock_decomposition(64,:,:),[3 2 1])*100  % DSD
Decomp_YSY=permute(x.oo_.shock_decomposition(64,:,:),[3 2 1])*100  % DSD



Decomp_q=permute(x.oo_.shock_decomposition(9,:,:),[3 2 1])*100  % q
Decomp_TOT=permute(x.oo_.shock_decomposition(10,:,:),[3 2 1])*100  % TOT
Decomp_mtotq=permute(x.oo_.shock_decomposition(69,:,:),[3 2 1])*100  % mtotq

Decomp_TShare=permute(x.oo_.shock_decomposition(54,:,:),[3 2 1])*100  % TShare
Decomp_XMY=permute(x.oo_.shock_decomposition(71,:,:),[3 2 1])*100  % XMY

