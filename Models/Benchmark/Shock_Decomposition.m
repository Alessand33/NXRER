
clear all;
close all;
clc;

x=load('dynare_Benchmark_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)

Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 
Decomp_DSD=permute(x.oo_.shock_decomposition(64,:,:),[3 2 1])*100  % DSD
Decomp_YSY=permute(x.oo_.shock_decomposition(61,:,:),[3 2 1])*100  % DSD


Decomp_q=permute(x.oo_.shock_decomposition(9,:,:),[3 2 1])*100  % q
Decomp_TOT=permute(x.oo_.shock_decomposition(10,:,:),[3 2 1])*100  % TOT
Decomp_mtotq=permute(x.oo_.shock_decomposition(69,:,:),[3 2 1])*100  % mtotq

Decomp_TShare=permute(x.oo_.shock_decomposition(54,:,:),[3 2 1])*100  % TShare
Decomp_XMY=permute(x.oo_.shock_decomposition(71,:,:),[3 2 1])*100  % XMY

Decomp_EXR=permute(x.oo_.shock_decomposition(41,:,:),[3 2 1])  % EXR
Decomp_IMR=permute(x.oo_.shock_decomposition(42,:,:),[3 2 1])  % IMR
Decomp_YR=permute(x.oo_.shock_decomposition(1,:,:),[3 2 1])  % IMR


%%  Parameter Values 
bet  = 0.99;          % time preference: real interest rate = 1 - 1/beta 
sig  = 6;             % CRRA coefficient 
th   = 4;             % Elasticity of demand across variety
rho  = 3.3;           % Elasticity of Substitution between Home and Foreign

rhob = 0.95;  % beta persistence 
fi   = 0.01;  % beta shock from C
sdb  = 0.00;  % volatility of shock to beta

xib  = 0.0001; % Bond cost parameter

rhozc = 0.98;  % Aggregate Productivity persistence
sdzc  = 0.0122; % SD(Z)

rhozd = 0.99;
sdzd = 0.0109;

rhoxic = 0.99;  % tariff persistence
sdci  = 0.0056; % SD(exi)

rhoxid = 0.99;  % tariff persistence
sddi  = 0.0476; % SD(exi)

rhoxia = 0.95;
sdtot = 0.01;
atotq = 0; % mean correction for totq
aysy = 0;  % mean correction for ysy

zeta = 0.5; % Shock to theta with RER

% No capital 
alp  = 0;             % capital share 
del  = 1;             % capital depreciation rate 

n1   = 0.025;        % stopper rate
N    = 0.20;            % # of exporters 
n0   = N/(1-N)*n1;  % entry rate  

%  modified parameters 
m    = th-1;  

sdeta  = 0.15;           % SD of firm specific shocks 

eta0 = norminv(1-n0,0,sdeta);  % marginal entrance productivity 
eta1 = norminv(n1,0,sdeta);  %  marginal exit productivity 

PsiT = exp(sdeta^2*m^2/2);
Psi0 = PsiT*(1-normcdf(eta0,m*sdeta^2,sdeta));
Psi1 = PsiT*(1-normcdf(eta1,m*sdeta^2,sdeta));
PsiX = N*Psi1 + (1-N)*Psi0;

Tr   = 0.10;             % Steady state IM/Y ratio 
L    = 1/4;              % Steady state Labor supply 

%% ================================================================================
%   Steady State Computation 
%  ================================================================================

xpar = [bet; sig; th; rho; sdeta; PsiT; Psi0; Psi1; PsiX; n1; ...
         n0; N; eta0; eta1; Tr; L; sdb];  % parameters used in SteadyS 
   
%---- Initial Values indm=1
EV0 = 1;
EV1 = 2;
dEV = 1;
a2 = 0.5; 
tau0 = 1;
tau1 = 0.1;

C    = 0.1; 
Ph   = 0.5; 
Pf   = 0.5; 
gam  = 0.35;

x=load('dynare_NoK02_results.mat');
NumberOfParameters = x.M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = deblank(x.M_.param_names(ii,:));
    eval([ paramname ' = x.M_.params(' int2str(ii) ');']);
end

% --- Steady State Computation 

x0= log([C;tau0;tau1;a2;gam;dEV;Ph;Pf]);
x = fsolve(@SS_NoK02,x0,[],xpar);

C    = exp(x(1)); 
tau0 = exp(x(2)); 
tau1 = exp(x(3));
a2   = exp(x(4)); 
gam  = exp(x(5));
dEV  = exp(x(6));
Ph   = exp(x(7));
Pf   = exp(x(8));

W = (1-gam)/gam*C/(1-L);
Lp  = L - (1-N)*n0*tau0 - N*(1-n1)*tau1;
EXY =  a2*Pf^(1-rho);
EV0 = ( 1/th*a2*(th*W/(th-1))^(1-th)*Psi0*Pf^(th-rho)*C - n0*W*tau0 + bet*n0*(dEV) )/(1-bet);
EV1 = dEV+EV0;

%% Dynare Simulation


q   = 1;
Z   = 1;
V   = bet;
B   = 0;
NXY = 0;
YN  = C;
YR  = C/Ph;
EXN = EXY*C;
IMN = EXN;
Px  = Pf*N^(1/(th-1));
Pm  = Pf*N^(1/(th-1));
EXR = EXN/Px;
IMR = IMN/Pm;
xi    = 1;  % iceberg cost normalized to 1
exi  = 1;
ez   = 1;
TOT = Px/Pm;
EXIMR = EXR/IMR;



Css = C;
Yss = log(YR);

Wh = W;
Wf = W;
EV0h = EV0;
EV0f = EV0;
EV1h = EV1;
EV1f = EV1;
eta0h = eta0;
eta0f = eta0;
eta1h = eta1;
eta1f = eta1;
Nh  = N;
Nf  = N;
Lh  = L;
Lf  = L;
Lph  = Lp;
Lpf = Lp;
Ph = Ph;
Pf = Pf;
Phs = Pf;
Pfs = Ph;
PsiXh = PsiX;
PsiXf = PsiX;
Ch = C;
Cf = C;
beth = bet;
betf = bet;
V  = V;
B = 0; 
q = 1;
NXY = 0;
YNh = YN;
YNf = YN;
YRh = YR;
YRf = YR;
EXN = EXN;
IMN = IMN;
Px = Px;
Pm = Pm;
EXR = EXR;
IMR = IMR; 
n0h = n0;
n0f = n0;
n1h = n1;
n1f = n1;
xih = 1;
xif = 1;
Zh = 1;
Zf = 1;
thh = th; 
thf = th;
xic = 1;
xid = 1;
TOT = 1;
EXIMR = 1;
TOTa = 1;
tshare = (EXN+IMN)/YN;
XMY  = (EXR+IMR)/YR;
NXYR = 1;
WEDGEDD = EXIMR/(Cf/Ch)*(TOTa/q)^rho;
WEDGEYY = EXIMR/(YRf/YRh)*(TOTa/q)^rho;
MD = IMR/YRh; 
YSY = YRf/YRh;
totq = TOTa/q;
ipdetrend = YRh/exp(Yss);
DSD = Cf/Ch;
Zc = 1;
Zd = 1;
By = 1;
EXIMN =EXN/IMN;
mtotq = exp(atotq);
mysy  = exp(aysy);
DSD_ADV = Cf/Ch;
mDSD_ADV = exp(aysy);
xia   = 1;
mtot = TOTa;



xpar=[bet sig th rho rhob fi xib rhozc sdzc rhoxic ...
      sdci zeta sdeta PsiT tau0 tau1 a2 gam Css rhoxid ...
      sddi sdzd rhozd Yss Psi0 Psi1 PsiX n1 n0 N ... 
      eta0 eta1 Tr L sdb atotq aysy rhoxia sdtot];

xss=[ Wh Wf EV0h EV0f EV1h EV1f exp(eta0h) exp(eta0f) exp(eta1h) exp(eta1f) ...
      Nh Nf Lh Lf Lph Lpf Ph Pfs Phs Pf ...
      PsiXh PsiXf Ch Cf beth betf V exp(B) q exp(NXY) ...
      YNh YNf YRh YRf EXN IMN Px Pm EXR IMR ... 
      n0h n0f n1h n1f xih xif Zh Zf thh thf ...
      xic xid TOT EXIMR TOTa tshare NXYR WEDGEDD WEDGEYY MD...
      YSY totq ipdetrend DSD_ADV Zc Zd By EXIMN mtotq xia ...
      XMY mDSD_ADV mysy mtot];

xss = log(xss);

save xpar;
save xss;

dynare dynare_NoK02_decomp;


x=load('dynare_NoK02_results.mat');

% % Loading shocks
ex_=[];
for shock_iter=1:x.M_.exo_nbr
    ex_=[ex_ x.oo_.SmoothedShocks.(deblank(x.M_.exo_names(shock_iter,:)))];
end
% save ex_;
% % use shocks only starting at t=2 due to t=1 being initial condition
ex1_ = ex_(1:end,:);
% 
% finding initial values of variables
y0=[];
for endo_iter=1:x.M_.endo_nbr
	y0 = [y0;
         x.oo_.SmoothedVariables.(deblank(x.M_.endo_names(endo_iter,:)))(1)];
end;
% save y0;
% 
ysteady_state = x.oo_.steady_state;
% save ysteady_state;

initial_condition_states = repmat(x.oo_.dr.ys,1,x.M_.maximum_lag);


% Variable printing to be used in Start_NoKcdshocks_feed.m 
for i=1:x.M_.endo_nbr
   name = x.M_.endo_names(i,:);
   name2 = '=ysteady_state(';
   name3 = ');';
   vnumb = i;
   fprintf('%s %s %d %s \n ',name,name2,vnumb,name3)
end;
disp(x.M_.endo_names(1:x.M_.endo_nbr,:));

% Parameters printing to be used in Start_NoKcdshocks_feed.m 
for i=1:x.M_.param_nbr
   name = x.M_.param_names(i,:);
   name2 = '=xpar1(';
   name3 = ');';
   vnumb = i;
   fprintf('%s %s %d %s \n ',name,name2,vnumb,name3)
end;

% Save Parameter values;
xpar1 = x.M_.params;
%save xpar1;

% use shocks only starting at t=2 due to t=1 being initial condition
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6
% Prod [1 2]  Trade [3 4]



ex_all = ex_(2:end,:);
ex_trade = ex_all;
ex_trade(:,[1 2 5]) = zeros(size(ex_all,1),3);
ex_prod = ex_all;
ex_prod(:,[3 4 5 6]) = zeros(size(ex_all,1),4);
ex_discount = ex_all;
ex_discount(:,[1 2 3 4 6]) = zeros(size(ex_all,1),5);
 
ex_noT = ex_all;
ex_noT(:,[3 4]) = zeros(size(ex_all,1),2);

ex_noTC = ex_all;
ex_noTC(:,[3]) = zeros(size(ex_all,1),1);

ex_noTD = ex_all;
ex_noTD(:,[4]) = zeros(size(ex_all,1),1);

dr = x.oo_.dr;


% y0 = initial_condition_states;
% 1. All Shocks + Initial
iorder = 1;
ex_tmp = ex_all;
yd_tmp = y0; 
% yd_tmp = ysteady_state;
Decomp2_A=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_A=Decomp2_A';
 
% 4. Trade Shock+Initial Trade shock
iorder = 1;
ex_tmp = ex_trade;
yd_tmp = y0; 
yd_tmp = ysteady_state;
yd_tmp([47 48 59 60],:) = y0([47 48 59 60],:); % xih xif xic xid 
Decomp2_T=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_T=Decomp2_T';


% 3. Productivity Shock+Initial
iorder = 1;
ex_tmp = ex_prod;
yd_tmp = y0; 
yd_tmp = ysteady_state;
yd_tmp([49 50 65 66],:) = y0([49 50 65 66],:); % zh zf zc zd 

Decomp2_P=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_P=Decomp2_P';

% 3. Discount Shock+Initial
iorder = 1;
ex_tmp = ex_discount;
yd_tmp = y0; 
yd_tmp = ysteady_state;
yd_tmp([15 16],:) = y0([15 16],:); % beth betf 
Decomp2_D=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_D=Decomp2_D';


% 3. Initial
iorder = 1;
ex_tmp = ex_discount*0;
yd_tmp = y0;
% yd_tmp = ysteady_state;
yd_tmp([47 48 59 60 49 50 65 66 15 16],:) = ysteady_state([47 48 59 60 49 50 65 66 15 16],:); % 
Decomp2_I=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_I=Decomp2_I';



% 3. No Trade shocks
iorder = 1;
ex_tmp = ex_noT;
yd_tmp = y0;
% yd_tmp = ysteady_state;
Decomp2_NoT=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_NoT=Decomp2_NoT';

% 3. No Common Trade shocks
iorder = 1;
ex_tmp = ex_noTC;
yd_tmp = y0;
% yd_tmp = ysteady_state;
Decomp2_NoTC=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_NoTC=Decomp2_NoTC';

% 3. No Diff Trade shocks
iorder = 1;
ex_tmp = ex_noTD;
yd_tmp = y0;
% yd_tmp = ysteady_state;
Decomp2_NoTD=simult_(yd_tmp,dr,ex_tmp,iorder);
Decomp2_NoTD=Decomp2_NoTD';



