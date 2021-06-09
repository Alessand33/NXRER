%% TABLE 5
% This program uses the benchmark model 
% It generates Excel files for simulated data:
% IR_SMPLE_All.xlsx: All shocks
% IR_SMPLE_NoB.xlsx: No beta shocks
% IR_SMPLE_NoXi.xlsx: No trade shocks
% IR_SMPLE_NoZ.xlsx: No productivity shocks
%
% Then using Stata, the Armington elasticities are estimated with these data set


clear all;
close all;
clc;
%%  Parameter Values 
bet  = 0.99;          % time preference: real interest rate = 1 - 1/beta 
sig  = 8;             % CRRA coefficient 
th   = 4;             % Elasticity of demand across variety
rho  = 3;           % Elasticity of Substitution between Home and Foreign

rhob = 0.95;  % beta persistence 
fi   = 0.002;  % beta shock from C
sdb  = 0.00;  % volatility of shock to beta

xib  = 0.0001; % Bond cost parameter

rhozc = 0.98;  % Aggregate Productivity persistence
sdzc  = 0.0122; % SD(Z)

rhozd = 0.99;
sdzd = 0.0109;

rhoxic = 0.9999;  % tariff persistence
sdci  = rho*0.01; % SD(exi)

rhoxid = 0.975;  % tariff persistence
sddi  = 0.0476; % SD(exi)

rhoxia = 0.95;
sdtot = 0.0225;
atotq = -0.1095; % mean correction for totq
aysy = -0.3326;  % mean correction for ysy

zeta = 0.56; % Shock to theta with RER

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

x=load('dynare_Benchmark_results.mat');
NumberOfParameters = x.M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = deblank(x.M_.param_names(ii,:));
    eval([ paramname ' = x.M_.params(' int2str(ii) ');']);
end

% Loading shocks
ex_=[];
for shock_iter=1:x.M_.exo_nbr
    ex_=[ex_ x.oo_.SmoothedShocks.(deblank(x.M_.exo_names(shock_iter,:)))];
end
save ex_;
Cov_ex = cov(ex_);
Var_ex = sqrt(diag(cov(ex_)));
Var_ex = Var_ex*Var_ex';
%Corr_ex = Cov_ex./Var_ex;
Corr_ex = Cov_ex;
CholCex = chol(Corr_ex); % Upper triangular  CholCex'*CholCex = Corr_ex 

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

% --- Steady State Computation 

x0= log([C;tau0;tau1;a2;gam;dEV;Ph;Pf]);
x = fsolve(@SS_Benchmark,x0,[],xpar);

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

dynare dynare_Benchmark_IR;


ysteady_state = oo_.steady_state;
% save ysteady_state;

ex_ = zeros(200,6);
% No Shocks
ex_all = ex_(1:end,:)*0;
% ex_all(2,2)=-0.1/sdzd;
% for i=3:size(ex_all,1);
%     ex_all(i,2)=ex_all(2,2)-rhozd*ex_all(2,2);
% end;
iorder = 1;
ex_tmp = ex_all;
yd_tmp = ysteady_state; 

% Variable printing to be used in Start_NoKcdshocks_feed.m 
for i=1:M_.endo_nbr
   name = M_.endo_names(i,:);
   name2 = '=ysteady_state(';
   name3 = ');';
   vnumb = i;
   fprintf('%s %s %d %s \n ',name,name2,vnumb,name3)
end;
disp(M_.endo_names(1:M_.endo_nbr,:));

% Parameters printing to be used in Start_NoKcdshocks_feed.m 
for i=1:M_.param_nbr
   name = M_.param_names(i,:);
   name2 = '=xpar1(';
   name3 = ');';
   vnumb = i;
   fprintf('%s %s %d %s \n ',name,name2,vnumb,name3)
end;

% Save Parameter values;
xpar1 = M_.params;
%save xpar1;

% use shocks only starting at t=2 due to t=1 being initial condition
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6
% Prod [1 2]  Trade [3 4]

dr = oo_.dr;


stream = RandStream.create('mt19937ar','seed',100);
RandStream.setGlobalStream(stream);

CholCex = eye(size(CholCex));  % no correlations 

%% All Shocks excluding TOT
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6

DropT = 100;
T = 100000+DropT;
N = 1;
Iter=1;

ex_tmp0 = randn(T*Iter,6*N);
for i=1:N;
    ex_tmp0(:,[(1+(i-1)*6):(i*6)]) = ex_tmp0(:,[(1+(i-1)*6):(i*6)])*CholCex; 
end;

IR_SMPLS = [];
CorXMAll=[];
CorXMDAll=[];
for i=1:N; 
    IR_SMPLE2=[];
    for j=1:Iter;
        ex_tmp = ex_tmp0([(1+(j-1)*T):(j*T)],[(1+(i-1)*6):(i*6)]);
        ex_tmp(:,6)=ex_tmp(:,6)*0; 
        IR_Simul=simult_(yd_tmp,dr,ex_tmp,iorder);
        IR_Simul=IR_Simul';
        IR_SMPLE = IR_Simul(DropT+2:end,[11 62 64 9 73]); % EXIMR TOTQ DSD Q YSY 
        IR_SMPLE2=[IR_SMPLE2;IR_SMPLE];
        [XCF1,lags1,bounds1]=crosscorr(IR_SMPLE(:,4),IR_SMPLE(:,1),12);
        [XCF2,lags2,bounds2]=crosscorr(IR_SMPLE(:,4),(IR_SMPLE(:,1)-IR_SMPLE(:,3)),12);
        CorXMAll=[CorXMAll;XCF1'];
        CorXMDAll=[CorXMDAll;XCF2'];
    end;
    IR_SMPLS=[IR_SMPLS,IR_SMPLE2];
end;


% ---- Writing in Excel ------
pnames = {'EXIMR1' 'TOTQ1' 'DSD1' 'Q1' 'YSY1'};
excelclear = nan(100001, 5);
xlswrite('IR_SMPLE_All.xlsx',excelclear,1,'A1:E100001');
xlswrite('IR_SMPLE_All.xlsx',pnames,1,'A1');
xlswrite('IR_SMPLE_All.xlsx',IR_SMPLS(1:100000,:),1,'A2');



%% No PROD SHOCKS: Excluding TOT Zc Zd
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6

IR_SMPLS = [];
CorXMNoZ=[];
CorXMDNoZ=[];
for i=1:N; 
    IR_SMPLE2=[];
    for j=1:Iter;
        ex_tmp = ex_tmp0([(1+(j-1)*T):(j*T)],[(1+(i-1)*6):(i*6)]);
        ex_tmp(:,6)=ex_tmp(:,6)*0;
        ex_tmp(:,1)=ex_tmp(:,1)*0;
        ex_tmp(:,2)=ex_tmp(:,2)*0;
        IR_Simul=simult_(yd_tmp,dr,ex_tmp,iorder);
        IR_Simul=IR_Simul';
        IR_SMPLE = IR_Simul(DropT+2:end,[11 62 64 9 73]); % EXIMR TOTQ DSD Q YSY 
        IR_SMPLE2=[IR_SMPLE2;IR_SMPLE];
        [XCF1,lags1,bounds1]=crosscorr(IR_SMPLE(:,4),IR_SMPLE(:,1),12);
        [XCF2,lags2,bounds2]=crosscorr(IR_SMPLE(:,4),(IR_SMPLE(:,1)-IR_SMPLE(:,3)),12);
        CorXMNoZ=[CorXMNoZ;XCF1'];
        CorXMDNoZ=[CorXMDNoZ;XCF2'];
    end;
    IR_SMPLS=[IR_SMPLS,IR_SMPLE2];
end;


% ---- Writing in Excel ------
pnames = {'EXIMR1' 'TOTQ1' 'DSD1' 'Q1' 'YSY1'};
excelclear = nan(100001, 5);
xlswrite('IR_SMPLE_NoZ.xlsx',excelclear,1,'A1:E100001');
xlswrite('IR_SMPLE_NoZ.xlsx',pnames,1,'A1');
xlswrite('IR_SMPLE_NoZ.xlsx',IR_SMPLS(1:100000,:),1,'A2');


%% No Trade Shocks: Excluding TOT Xic Xid
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6

IR_SMPLS = [];
CorXMNoXi=[];
CorXMDNoXi=[];
for i=1:N; 
    IR_SMPLE2=[];
    for j=1:Iter;
        ex_tmp = ex_tmp0([(1+(j-1)*T):(j*T)],[(1+(i-1)*6):(i*6)]);
        ex_tmp(:,6)=ex_tmp(:,6)*0;
        ex_tmp(:,3)=ex_tmp(:,3)*0;
        ex_tmp(:,4)=ex_tmp(:,4)*0;
        IR_Simul=simult_(yd_tmp,dr,ex_tmp,iorder);
        IR_Simul=IR_Simul';
        IR_SMPLE = IR_Simul(DropT+2:end,[11 62 64 9 73]); % EXIMR TOTQ DSD Q YSY 
        IR_SMPLE2=[IR_SMPLE2;IR_SMPLE];
        [XCF1,lags1,bounds1]=crosscorr(IR_SMPLE(:,4),IR_SMPLE(:,1),12);
        [XCF2,lags2,bounds2]=crosscorr(IR_SMPLE(:,4),(IR_SMPLE(:,1)-IR_SMPLE(:,3)),12);
        CorXMNoXi=[CorXMNoXi;XCF1'];
        CorXMDNoXi=[CorXMDNoXi;XCF2'];
    end;
    IR_SMPLS=[IR_SMPLS,IR_SMPLE2];
end;


% ---- Writing in Excel ------
pnames = {'EXIMR1' 'TOTQ1' 'DSD1' 'Q1' 'YSY1'};
excelclear = nan(100001, 5);
xlswrite('IR_SMPLE_NoXi.xlsx',excelclear,1,'A1:E100001');
xlswrite('IR_SMPLE_NoXi.xlsx',pnames,1,'A1');
xlswrite('IR_SMPLE_NoXi.xlsx',IR_SMPLS(1:100000,:),1,'A2');


%% No Beta Shocks: Excluding TOT Beta
% ezc ezd exic exid ebf etot
% 1    2   3    4   5  6

IR_SMPLS = [];
CorXMNoB=[];
CorXMDNoB=[];
for i=1:N; 
    IR_SMPLE2=[];
    for j=1:Iter;
        ex_tmp = ex_tmp0([(1+(j-1)*T):(j*T)],[(1+(i-1)*6):(i*6)]);
        ex_tmp(:,6)=ex_tmp(:,6)*0;
        ex_tmp(:,5)=ex_tmp(:,5)*0;
        
        IR_Simul=simult_(yd_tmp,dr,ex_tmp,iorder);
        IR_Simul=IR_Simul';
        IR_SMPLE = IR_Simul(DropT+2:end,[11 62 64 9 73]); % EXIMR TOTQ DSD Q YSY 
        IR_SMPLE2=[IR_SMPLE2;IR_SMPLE];
        [XCF1,lags1,bounds1]=crosscorr(IR_SMPLE(:,4),IR_SMPLE(:,1),12);
        [XCF2,lags2,bounds2]=crosscorr(IR_SMPLE(:,4),(IR_SMPLE(:,1)-IR_SMPLE(:,3)),12);
        CorXMNoB=[CorXMNoB;XCF1'];
        CorXMDNoB=[CorXMDNoB;XCF2'];
    end;
    IR_SMPLS=[IR_SMPLS,IR_SMPLE2];
end;


% ---- Writing in Excel ------
pnames = {'EXIMR1' 'TOTQ1' 'DSD1' 'Q1' 'YSY1'};
excelclear = nan(100001, 5);
xlswrite('IR_SMPLE_NoB.xlsx',excelclear,1,'A1:E100001');
xlswrite('IR_SMPLE_NoB.xlsx',pnames,1,'A1');
xlswrite('IR_SMPLE_NoB.xlsx',IR_SMPLS(1:100000,:),1,'A2');



