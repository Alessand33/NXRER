clc;
close all;
%%-----------------------------------------------------------------------------------------------------------------------------------------------
var YRh, YRf, Ch, Cf, Lh, Lf, Lph, Lpf, q, TOT, 
    EXIMR, NXY, Nh, Nf, beth, betf, Wh, Wf, EV0h, EV0f, 
    EV1h, EV1f, eta0h, eta0f, eta1h, eta1f, Ph, Pfs, Phs, Pf, 
    PsiXh, PsiXf, V, B, YNh, YNf, EXN, IMN, Px, Pm, 
    EXR, IMR, n0h, n0f, n1h, n1f, xih, xif, Zh, Zf, 
    thh, thf, TOTa, tshare, NXYR, WEDGEDD, WEDGEYY, MD, xic, xid, 
    YSY, totq, ipdetrend, DSD_ADV, Zc, Zd, By, EXIMN, mtotq, xia, 
    XMY, mDSD_ADV, mysy, mtot;

varexo ezc, ezd, exic, exid, ebf, etot;

load xss;
load xpar;
%%-----------------------------------------------------------------------------------------------------------------------------------------------
parameters bet, sig, th, rho, rhob, fi, xib, rhozc, sdzc, rhoxic, 
           sdci, zeta, sdeta, PsiT, tau0, tau1, a2, gam, Css, rhoxid, 
           sddi, sdzd, rhozd, Yss, Psi0, Psi1, PsiX, n1, n0, N, 
           eta0, eta1, Tr, L, sdb, atotq, aysy, rhoxia, sdtot;

bet   = xpar(1);
sig   = xpar(2);
th    = xpar(3);
rho   = xpar(4);
rhob  = xpar(5);
fi    = xpar(6);
xib   = xpar(7);
rhozc = xpar(8);
sdzc  = xpar(9);
rhoxic = xpar(10);
sdci  = xpar(11);
zeta  = xpar(12);
sdeta = xpar(13);
PsiT  = xpar(14);
tau0  = xpar(15);
tau1  = xpar(16);
a2    = xpar(17);
gam   = xpar(18);
Css   = xpar(19);
rhoxid = xpar(20);
sddi   = xpar(21);
sdzd = xpar(22);
rhozd = xpar(23);
Yss = xpar(24);
Psi0  = xpar(25); 
Psi1 = xpar(26);
PsiX = xpar(27);
n1 = xpar(28);
n0 = xpar(29); 
N  = xpar(30); 
eta0 = xpar(31); 
eta1 = xpar(32); 
Tr = xpar(33); 
L = xpar(34);
sdb = xpar(35);
atotq = xpar(36);
aysy = xpar(37);
rhoxia = xpar(38);
sdtot = xpar(39);
%%-----------------------------------------------------------------------------------------------------------------------------------------------
model; 

% 1 Wh
(1-gam)/gam*exp(Ch-Wh) = ( 1 - exp(Lh) );   
% 2 Wf
(1-gam)/gam*exp(Cf-Wf) = ( 1 - exp(Lf) );
% 3 EV0h 
exp(EV0h) = 1/exp(thf)*a2*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta0h,(exp(thf)-1)*sdeta^2,sdeta))*exp( exp(thf)*q+(exp(thf)-1)*Zh + (exp(thf)-rho)*Phs + Cf)
            - exp(n0h+Wh)*tau0 + exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )*( exp(EV0h(+1)) + exp(n0h)*( exp(EV1h(+1))-exp(EV0h(+1)) ) );
% 4 EV0f
exp(EV0f) = 1/exp(thh)*a2*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta0f,(exp(thh)-1)*sdeta^2,sdeta))*exp(-exp(thh)*q+(exp(thh)-1)*Zf + (exp(thh)-rho)*Pf + Ch) 
            - exp(n0f+Wf)*tau0 + exp(betf)*exp( sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )*( exp(EV0f(+1)) + exp(n0f)*( exp(EV1f(+1))-exp(EV0f(+1)) ) );
% 5 EV1h
exp(EV1h) = 1/exp(thf)*a2*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta1h,(exp(thf)-1)*sdeta^2,sdeta))*exp( exp(thf)*q+(exp(thf)-1)*Zh + (exp(thf)-rho)*Phs + Cf)
            - (1-exp(n1h))*exp(Wh)*tau1 + exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )*( exp(EV0h(+1)) + (1-exp(n1h))*( exp(EV1h(+1))-exp(EV0h(+1)) ) );
% 6 EV1f
exp(EV1f) = 1/exp(thh)*a2*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta1f,(exp(thh)-1)*sdeta^2,sdeta))*exp(-exp(thh)*q+(exp(thh)-1)*Zf + (exp(thh)-rho)*Pf + Ch) 
            - (1-exp(n1f))*exp(Wf)*tau1 + exp(betf)*exp( sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )*( exp(EV0f(+1)) + (1-exp(n1f))*( exp(EV1f(+1))-exp(EV0f(+1)) ) );
% 7 eta0h
exp(Wh)*tau0 = 1/exp(thf)*a2*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp( exp(thf)*q+(exp(thf)-1)*(Zh+eta0h) + (exp(thf)-rho)*Phs + Cf) 
                + exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )*( exp(EV1h(+1))-exp(EV0h(+1)) );
% 8 eta0f
exp(Wf)*tau0 = 1/exp(thh)*a2*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(-exp(thh)*q+(exp(thh)-1)*(Zf+eta0f) + (exp(thh)-rho)*Pf  + Ch) 
                + exp(betf)*exp( sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )*( exp(EV1f(+1))-exp(EV0f(+1)) );
% 9 eta1h
exp(Wh)*tau1 = 1/exp(thf)*a2*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp( exp(thf)*q+(exp(thf)-1)*(Zh+eta1h) + (exp(thf)-rho)*Phs + Cf) 
                + exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )*( exp(EV1h(+1))-exp(EV0h(+1)) );
% 10 eta0f
exp(Wf)*tau1 = 1/exp(thh)*a2*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(-exp(thh)*q+(exp(thh)-1)*(Zf+eta1f) + (exp(thh)-rho)*Pf  + Ch) 
                + exp(betf)*exp( sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )*( exp(EV1f(+1))-exp(EV0f(+1)) );
% 11 Nh
exp(Nh) = (1-normcdf(eta0h,0,sdeta))*(1-exp(Nh(-1)))+(1-normcdf(eta1h,0,sdeta))*exp(Nh(-1));
% 12 Nf 
exp(Nf) = (1-normcdf(eta0f,0,sdeta))*(1-exp(Nf(-1)))+(1-normcdf(eta1f,0,sdeta))*exp(Nf(-1));
% 13 Lh
exp(Lh) = exp(Lph) + (1-exp(Nh(-1)))*(1-normcdf(eta0h,0,sdeta))*tau0 + exp(Nh(-1))*(1-normcdf(eta1h,0,sdeta))*tau1;
% 14 Lf
exp(Lf) = exp(Lpf) + (1-exp(Nf(-1)))*(1-normcdf(eta0f,0,sdeta))*tau0 + exp(Nf(-1))*(1-normcdf(eta1f,0,sdeta))*tau1;
% 15 Lph
exp(Lph) = (th*exp(Wh)/(th-1))^(-1)*exp((1-rho)*Ph+Ch) + a2*(exp(thf)*exp(Wh)/(exp(thf)-1))^(-1)*exp( q+(1-rho)*Phs+Cf);
% 16 Lpf
exp(Lpf) = (th*exp(Wf)/(th-1))^(-1)*exp((1-rho)*Pfs+Cf) + a2*(exp(thh)*exp(Wf)/(exp(thh)-1))^(-1)*exp(-q+(1-rho)*Pf+Ch);
% 17 Ph 
exp((1-th)*Ph)  = (th*exp(Wh-Zh)/(th-1))^(1-th)*PsiT;
% 18 Pfs 
exp((1-th)*Pfs) = (th*exp(Wf-Zf)/(th-1))^(1-th)*PsiT;
% 19 Phs 
exp((1-exp(thf))*(Phs+q)) = exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh-Zh)/(exp(thf)-1))^(1-exp(thf))*exp(PsiXh);
% 20 Pf 
exp((1-exp(thh))*(Pf-q)) = exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf-Zf)/(exp(thh)-1))^(1-exp(thh))*exp(PsiXf);
% 21 PsiXh
exp(PsiXh) = exp(Nh(-1))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta1h,(exp(thf)-1)*sdeta^2,sdeta)) + (1-exp(Nh(-1)))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta0h,(exp(thf)-1)*sdeta^2,sdeta));
% 22 PsiXf
exp(PsiXf) = exp(Nf(-1))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta1f,(exp(thh)-1)*sdeta^2,sdeta)) + (1-exp(Nf(-1)))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta0f,(exp(thh)-1)*sdeta^2,sdeta));
% 23 Ch 
1 = exp((1-rho)*Ph) + a2*exp((1-rho)*Pf);
% 24 Cf
1 = exp((1-rho)*Pfs) + a2*exp((1-rho)*Phs);
% 25 beth
beth = (1-rhob)*log(bet)+fi*log(Css) + rhob*beth(-1)-fi*Ch - sdb*ebf/2;
% 26 betf
betf = (1-rhob)*log(bet)+fi*log(Css) + rhob*betf(-1)-fi*Cf + sdb*ebf/2;
% 27 V
exp(V)*(1+xib*B*exp(V-YNh)) = exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) );
% 28 B 
exp(Ch) + exp(V)*B = exp(Wh+Lph) + B(-1) + 1/th*exp((1-rho)*Ph+Ch) + 1/exp(thf)*a2*exp(q+(1-rho)*Phs+Cf);
% 29 q
exp(beth)*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )/(1+xib*B*exp(V-YNh))
    = exp(betf)*exp( q-q(+1) + sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )/(1-xib*B*exp(V-YNf-q));
% 30 NXY
NXY = (exp(EXN)-exp(IMN))/exp(YNh);
% 31 YNh
exp(YNh) = exp( (1-rho)*Ph+Ch) + a2*exp( q + (1-rho)*Phs + Cf );
% 32 YNf
exp(YNf) = exp( (1-rho)*Pfs+Cf) + a2*exp(-q + (1-rho)*Pf + Ch );
% 33 YRh 
exp(YRh) = exp(YNh-Ph);
% 34 YRf 
exp(YRf) = exp(YNf-Pfs);
% 35 EXN 
EXN = log(a2)+q+(1-rho)*Phs + Cf;
% 36 IMN 
IMN = log(a2)+(1-rho)*Pf + Ch;
% 37 Px 
exp(Px) = exp(Phs+q-xif+1/(exp(thf)-1)*Nh) + sdtot*etot/2;
% 38 Pm
exp(Pm) = exp(Pf-xih+1/(exp(thh)-1)*Nf) - sdtot*etot/2;
% 39 EXR
exp(EXR) = exp(EXN-Px);
% 40 IMR
exp(IMR) = exp(IMN-Pm);
% 41 n0h 
1-normcdf(eta0h,0,sdeta) = exp(n0h);
% 42 n0f
1-normcdf(eta0f,0,sdeta) = exp(n0f);
% 43 n1h
normcdf(eta1h,0,sdeta) = exp(n1h);
% 44 n1f
normcdf(eta1f,0,sdeta) = exp(n1f);
% 45 xih
xih = xic+0.5*xid + xia;
% 46 xif
xif = xic-0.5*xid;
% 47 Zh
Zc = rhozc*Zc(-1) + sdzc*ezc;
% 48 Zf
Zd = rhozd*Zd(-1) + sdzd*ezd;
% 49
Zh = Zc+0.5*Zd;
% 50
Zf = Zc-0.5*Zd;

% 51 thh
exp(thh) = th*exp(zeta*q);
% 52 thf
exp(thf) = th*exp(-zeta*q);
% 53 TOT
TOT = Pm-Px;
% 54 EXIMR
EXIMR = EXR-IMR;
% 55 TOTa
TOTa = Pm-Px;
% 56 tshare
tshare=(exp(EXN)+exp(IMN))/exp(YNh);
% 57 NXYR
NXYR = (exp(EXR)-exp(IMR))/exp(YRh);
% 58 WEDGEDD  
WEDGEDD = EXIMR - (Cf - Ch) - rho*totq;
% 59 WEDGEYY  
WEDGEYY = EXIMR - (YRf - YRh) - rho*totq;
% 60 IMPORTSHARE
MD = IMR - YRh;
% 61 xic
xic = rhoxic*xic(-1) + sdci*exic/rho;
% 62 xid
xid = rhoxid*xid(-1) + sddi*exid/rho;
% 63 YSY
YSY = YRf - YRh;
% 64 totq
totq = TOT+q;
% 65 ipdetrend
ipdetrend = YRh - Yss;
% 66 DSD_ADV
DSD_ADV = Cf- Ch;
% 67 Assets to Nominal GDP
By = B/exp(YNh);
% 68 Nominal EXPORT IMPORT RATIO
EXIMN = EXN-IMN;

% 69 mtotq
mtotq = totq + atotq;

%70 mDSD_ADV
mDSD_ADV = DSD_ADV + aysy;

%71 xia
xia = rhoxia*xia(-1);

%72 XMY 
exp(XMY) = ( exp(EXR) + exp(IMR) )/exp(YRh);

% 73 mysy
mysy = YSY+aysy;


% 74 mtot
mtot = TOT ;

end;
%---------------------------------------------
%  -------------  END OF MODEL   -------------
%---------------------------------------------

initval;

Wh    = xss(1);
Wf    = xss(2);
EV0h  = xss(3);
EV0f  = xss(4);
EV1h  = xss(5);
EV1f  = xss(6);
eta0h = xss(7);
eta0f = xss(8);
eta1h = xss(9);
eta1f = xss(10);
Nh    = xss(11);
Nf    = xss(12);
Lh    = xss(13);
Lf    = xss(14);
Lph   = xss(15);
Lpf   = xss(16);
Ph    = xss(17);
Pfs   = xss(18);
Phs   = xss(19);
Pf    = xss(20);
PsiXh = xss(21);
PsiXf = xss(22);
Ch    = xss(23);
Cf    = xss(24);
beth  = xss(25);
betf  = xss(26);
V     = xss(27);
B     = xss(28);
q     = xss(29);
NXY   = xss(30);
YNh   = xss(31);
YNf   = xss(32);
YRh   = xss(33);
YRf   = xss(34);
EXN   = xss(35);
IMN   = xss(36);
Px    = xss(37);
Pm    = xss(38);
EXR   = xss(39);
IMR   = xss(40);
n0h   = xss(41);
n0f   = xss(42);
n1h   = xss(43);
n1f   = xss(44);
xih   = xss(45);
xif   = xss(46);
Zh    = xss(47);
Zf    = xss(48);
thh   = xss(49);
thf   = xss(50);
xic  = xss(51); 
xid  = xss(52); 
TOT   = xss(53);
EXIMR = xss(54);
TOTa   = xss(55);
tshare = xss(56);
NXYR   = xss(57);
WEDGEDD = xss(58);
WEDGEYY = xss(59);
MD = xss(60);
YSY = xss(61);
totq = xss(62);
ipdetrend = xss(63);
DSD_ADV = xss(64);
Zc =xss(65);
Zd=xss(66);
By = xss(67);
EXIMN = xss(68);
mtotq = xss(69);
xia  = xss(70);
XMY  = xss(71);
mDSD_ADV = xss(72);
ebf=0;
mysy = xss(73);
mtot =xss(74);
end;


%-----------------------
%----    SHOCKS   ------
%-----------------------
shocks;

var ezc = 1;  % aggregate productivity 
var ezd = 1;
var ezc,ezd = 0.0;
var exic = 1;  % Marginal Trade Costs
var exid = 1;
var exic,exid = 0.0;
var ebf = 1;    % beta shocks
var etot = 1;

end;

%-----------------------
%----    SOLVE THE MODEL   ------
%-----------------------
steady;
check;
resid(1);

stoch_simul(order=1, graph) EXIMN EXIMR mtot;

estimated_params;

rhozc, 0.98, uniform_pdf, 0.9,1.00;
rhozd, 0.95, uniform_pdf, 0.9,1.00;
rhoxic, 0.98, uniform_pdf, 0.90,1;
rhoxid, 0.95, uniform_pdf, 0.90,1;
rhob, 0.95, uniform_pdf, , , 0.9,1;
sdci , 0.01, inv_gamma_pdf, 0.01, 0.2;
sddi , 0.01, inv_gamma_pdf, 0.01, 0.2;
sdzc , 0.01, inv_gamma_pdf, 0.01, 0.02;
sdzd , 0.007, inv_gamma_pdf, 0.01, 0.02;
sdb,  0.001, inv_gamma_pdf, 0.001, 0.02;
sdtot, 0.01, inv_gamma_pdf, 0.01, 0.05;
zeta, 0.5, normal_pdf, 0.5,0.1;
fi, 0.01, inv_gamma_pdf, 0.01, 0.1;
rho, 3.3, uniform_pdf, , , 1, 3.99;
sig, 6, uniform_pdf, , , 1, 8;
atotq, 0, normal_pdf, 0, 0.2;
aysy,  0, normal_pdf, 0, 0.2;


end;

estimated_params_init(use_calibration);
end;

varobs EXIMR XMY mysy ipdetrend mtotq mtot;
estimation(datafile=Data02_world, mode_compute = 6, nobs=139,mh_replic=2500,mh_nblocks=3,mh_jscale=0.3, plot_priors = 0, nodisplay);
   