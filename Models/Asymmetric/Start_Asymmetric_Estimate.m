%%   Asymmetric Countries
%    Foreign relative to home is set with productivity in Foreign, +zgap,
%    to get the target relative size, relF, Yfq/Yh in nominal values.
%    
%    all the other parameters are the same as the benchmark. (the steady
%    state values for the targets are slightly different but minimal, e.g., L, N, EXY.


clear all;
close all;
clc;
%%  Parameter Values 
bet  = 0.99;          % time preference: real interest rate = 1 - 1/beta 
sig  = 6;             % CRRA coefficient 
th   = 4;             % Elasticity of demand across variety
rho  = 3.3;           % Elasticity of Substitution between Home and Foreign

rhob = 0.95;  % beta persistence 
fi   = 0.00;  % beta shock from C
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


relF = 1.4; % Foreign size relative to home nominal YfQ/Yh
zgap = 0.0;   % relative productivity (to be corrected to match relF)
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

a2h = a2;
a2f = a2;
tau0h = tau0;
tau0f = tau0;
tau1h = tau1;
tau1f = tau1;
zgap = 0.0;
Cssh = Css;
Cssf = Css;

xpar=[bet sig th rho rhob fi xib rhozc sdzc rhoxic ...
      sdci zeta sdeta PsiT tau0h tau1h a2h gam Cssh rhoxid ...
      sddi sdzd rhozd Yss Psi0 Psi1 PsiX n1 n0 N ... 
      eta0 eta1 Tr L sdb atotq aysy rhoxia sdtot a2f ...
      tau0f tau1f relF Cssf];


xss=[ Wh Wf EV0h EV0f EV1h EV1f exp(eta0h) exp(eta0f) exp(eta1h) exp(eta1f) ...
      Nh Nf Lh Lf Lph Lpf Ph Pfs Phs Pf ...
      PsiXh PsiXf Ch Cf beth betf V exp(B) q exp(NXY) ...
      YNh YNf YRh YRf EXN IMN Px Pm EXR IMR ... 
      n0h n0f n1h n1f xih xif Zh Zf thh thf ...
      xic xid TOT EXIMR TOTa tshare exp(NXYR) WEDGEDD WEDGEYY MD...
      YSY totq ipdetrend DSD_ADV Zc Zd By EXIMN mtotq xia ...
      XMY mDSD_ADV mysy mtot exp(zgap)];

xss = log(xss);

save xpar;
save xss;


% --- Steady State Computation 

%x0= log([C;tau0;tau1;a2;gam;dEV;Ph;Pf]);
x = fsolve(@Steady_Asymmetric,xss,[],xpar);

Wh    = exp(x(1));
Wf    = exp(x(2));
EV0h  = exp(x(3));
EV0f  = exp(x(4));
EV1h  = exp(x(5));
EV1f  = exp(x(6));
eta0h = (x(7));
eta0f = (x(8));
eta1h = (x(9));
eta1f = (x(10));
Nh    = exp(x(11));
Nf    = exp(x(12));
Lh    = exp(x(13));
Lf    = exp(x(14));
Lph   = exp(x(15));
Lpf   = exp(x(16));
Ph    = exp(x(17));
Pfs   = exp(x(18));
Phs   = exp(x(19));
Pf    = exp(x(20));
PsiXh = exp(x(21));
PsiXf = exp(x(22));
Ch    = exp(x(23));
Cf    = exp(x(24));
beth  = exp(x(25));
betf  = exp(x(26));
V     = exp(x(27));
B     = (x(28));
q     = exp(x(29));
NXY   = (x(30));
YNh   = exp(x(31));
YNf   = exp(x(32));
YRh   = exp(x(33));
YRf   = exp(x(34));
EXN   = exp(x(35));
IMN   = exp(x(36));
Px    = exp(x(37));
Pm    = exp(x(38));
EXR   = exp(x(39));
IMR   = exp(x(40));
n0h   = exp(x(41));
n0f   = exp(x(42));
n1h   = exp(x(43));
n1f   = exp(x(44));
xih   = exp(x(45));
xif   = exp(x(46));
Zh    = exp(x(47));
Zf    = exp(x(48));
thh   = exp(x(49));
thf   = exp(x(50));
xic   = exp(x(51)); 
xid   = exp(x(52)); 
TOT   = exp(x(53));
EXIMR = exp(x(54));
TOTa   = exp(x(55));
tshare = (x(56));
NXYR   = (x(57));
WEDGEDD = exp(x(58));
WEDGEYY = exp(x(59));
MD = exp(x(60));
YSY = exp(x(61));
totq = exp(x(62));
ipdetrend = exp(x(63));
DSD_ADV = exp(x(64));
Zc = exp(x(65));
Zd = exp(x(66));
By = (x(67));
EXIMN = exp(x(68));
mtotq = exp(x(69));
xia  = exp(x(70));
XMY  = exp(x(71));
mDSD_ADV = exp(x(72));
mysy = exp(x(73));
mtot =exp(x(74));
zgap =(x(75));

Cssh = Ch;
Cssf = Cf;
Yssh = log(YRh);


xpar=[bet sig th rho rhob fi xib rhozc sdzc rhoxic ...
      sdci zeta sdeta PsiT tau0h tau1h a2h gam Cssh rhoxid ...
      sddi sdzd rhozd Yssh Psi0 Psi1 PsiX n1 n0 N ... 
      eta0 eta1 Tr L sdb atotq aysy rhoxia sdtot a2f ...
      tau0f tau1f zgap Cssf];


xss=[ Wh Wf EV0h EV0f EV1h EV1f exp(eta0h) exp(eta0f) exp(eta1h) exp(eta1f) ...
      Nh Nf Lh Lf Lph Lpf Ph Pfs Phs Pf ...
      PsiXh PsiXf Ch Cf beth betf V exp(B) q exp(NXY) ...
      YNh YNf YRh YRf EXN IMN Px Pm EXR IMR ... 
      n0h n0f n1h n1f xih xif Zh Zf thh thf ...
      xic xid TOT EXIMR TOTa exp(tshare) exp(NXYR) WEDGEDD WEDGEYY MD...
      YSY totq ipdetrend DSD_ADV Zc Zd exp(By) EXIMN mtotq xia ...
      XMY mDSD_ADV mysy mtot];

xss = log(xss);

save xpar;
save xss;

dynare dynare_Asymmetric;
