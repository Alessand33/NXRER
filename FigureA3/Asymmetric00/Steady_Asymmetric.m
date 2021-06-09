function F = Steady_Asymmetric(xssasym,xparasym);

% tau0h --> eta0h, tau1h--> eta1h, tau0f --> eta0f, tau1f --> eta1f
% a2 --> IMN-Ch 
% zgap --> relF
% gam --> L
% 

bet   = xparasym(1);
sig   = xparasym(2);
th    = xparasym(3);
rho   = xparasym(4);
rhob  = xparasym(5);
fi    = xparasym(6);
xib   = xparasym(7);
rhozc = xparasym(8);
sdzc  = xparasym(9);
rhoxic = xparasym(10);
sdci  = xparasym(11);
zeta  = xparasym(12);
sdeta = xparasym(13);
PsiT  = xparasym(14);
tau0h  = xparasym(15);
tau1h  = xparasym(16);
a2h    = xparasym(17);
gam   = xparasym(18);
Cssh   = xparasym(19);
rhoxid = xparasym(20);
sddi   = xparasym(21);
sdzd = xparasym(22);
rhozd = xparasym(23);
Yss = xparasym(24);
Psi0  = xparasym(25); 
Psi1 = xparasym(26);
PsiX = xparasym(27);
n1 = xparasym(28);
n0 = xparasym(29); 
N  = xparasym(30); 
eta0 = xparasym(31); 
eta1 = xparasym(32); 
Tr = xparasym(33); 
L = xparasym(34);
sdb = xparasym(35);
atotq = xparasym(36);
aysy = xparasym(37);
rhoxia = xparasym(38);
sdtot = xparasym(39);
a2f = xparasym(40);
tau0f = xparasym(41);
tau1f = xparasym(42);
relF  = xparasym(43);
Cssf  = xparasym(44);

ezc = 0;  % aggregate productivity 
ezd = 0;
exic = 0;  % Marginal Trade Costs
exid = 0;
ebf = 0;    % beta shocks
etot = 0;

Wh    = xssasym(1);
Wf    = xssasym(2);
EV0h  = xssasym(3);
EV0f  = xssasym(4);
EV1h  = xssasym(5);
EV1f  = xssasym(6);
eta0h = xssasym(7);
eta0f = xssasym(8);
eta1h = xssasym(9);
eta1f = xssasym(10);
Nh    = xssasym(11);
Nf    = xssasym(12);
Lh    = xssasym(13);
Lf    = xssasym(14);
Lph   = xssasym(15);
Lpf   = xssasym(16);
Ph    = xssasym(17);
Pfs   = xssasym(18);
Phs   = xssasym(19);
Pf    = xssasym(20);
PsiXh = xssasym(21);
PsiXf = xssasym(22);
Ch    = xssasym(23);
Cf    = xssasym(24);
beth  = xssasym(25);
betf  = xssasym(26);
V     = xssasym(27);
B     = xssasym(28);
q     = xssasym(29);
NXY   = xssasym(30);
YNh   = xssasym(31);
YNf   = xssasym(32);
YRh   = xssasym(33);
YRf   = xssasym(34);
EXN   = xssasym(35);
IMN   = xssasym(36);
Px    = xssasym(37);
Pm    = xssasym(38);
EXR   = xssasym(39);
IMR   = xssasym(40);
n0h   = xssasym(41);
n0f   = xssasym(42);
n1h   = xssasym(43);
n1f   = xssasym(44);
xih   = xssasym(45);
xif   = xssasym(46);
Zh    = xssasym(47);
Zf    = xssasym(48);
thh   = xssasym(49);
thf   = xssasym(50);
xic  = xssasym(51); 
xid  = xssasym(52); 
TOT   = xssasym(53);
EXIMR = xssasym(54);
TOTa   = xssasym(55);
tshare = xssasym(56);
NXYR   = xssasym(57);
WEDGEDD = xssasym(58);
WEDGEYY = xssasym(59);
MD = xssasym(60);
YSY = xssasym(61);
totq = xssasym(62);
ipdetrend = xssasym(63);
DSD_ADV = xssasym(64);
Zc =xssasym(65);
Zd=xssasym(66);
By = xssasym(67);
EXIMN = xssasym(68);
mtotq = xssasym(69);
xia  = xssasym(70);
XMY  = xssasym(71);
mDSD_ADV = xssasym(72);
ebf=0;
mysy = xssasym(73);
mtot =xssasym(74);
zgap = xssasym(75);
%zgap = 0.1;

% 1 Wh
F(1)  = (1-gam)/gam*exp(Ch-Wh) - ( 1 - exp(Lh) );
% 2 Wf
F(2)  = (1-gam)/gam*exp(Cf-Wf) - ( 1 - exp(Lf) );

% 3 EV0h 
F(3)  = 1/exp(thf)*a2f*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta0h,(exp(thf)-1)*sdeta^2,sdeta))*exp( exp(thf)*q+(exp(thf)-1)*Zh + (exp(thf)-rho)*Phs + Cf) ...
            - exp(n0h+Wh)*tau0h + exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) )*( exp(EV0h) + exp(n0h)*( exp(EV1h)-exp(EV0h) ) ) - exp(EV0h);
% 4 EV0f
F(4)  = 1/exp(thh)*a2h*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta0f,(exp(thh)-1)*sdeta^2,sdeta))*exp(-exp(thh)*q+(exp(thh)-1)*Zf + (exp(thh)-rho)*Pf + Ch) ...
            - exp(n0f+Wf)*tau0f + exp(betf)*exp( sig*(Cf-Cf) + (1-gam)*(1-sig)*(Wf-Wf) )*( exp(EV0f) + exp(n0f)*( exp(EV1f)-exp(EV0f) ) ) - exp(EV0f);

% 5 EV1h
F(5)  = 1/exp(thf)*a2f*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta1h,(exp(thf)-1)*sdeta^2,sdeta))*exp( exp(thf)*q+(exp(thf)-1)*Zh + (exp(thf)-rho)*Phs + Cf) ...
            - (1-exp(n1h))*exp(Wh)*tau1h + exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) )*( exp(EV0h) + (1-exp(n1h))*( exp(EV1h)-exp(EV0h) ) )  - exp(EV1h);
% 6 EV1f
F(6)  = 1/exp(thh)*a2h*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta1f,(exp(thh)-1)*sdeta^2,sdeta))*exp(-exp(thh)*q+(exp(thh)-1)*Zf + (exp(thh)-rho)*Pf + Ch) ...
            - (1-exp(n1f))*exp(Wf)*tau1f + exp(betf)*exp( sig*(Cf-Cf) + (1-gam)*(1-sig)*(Wf-Wf) )*( exp(EV0f) + (1-exp(n1f))*( exp(EV1f)-exp(EV0f) ) ) - exp(EV1f);
% 7 eta0h
F(7)  = 1/exp(thf)*a2f*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp( exp(thf)*q+(exp(thf)-1)*(Zh+eta0h) + (exp(thf)-rho)*Phs + Cf) ...
                + exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) )*( exp(EV1h)-exp(EV0h) ) - exp(Wh)*tau0h;
% 8 eta0f
F(8)  = 1/exp(thh)*a2h*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(-exp(thh)*q+(exp(thh)-1)*(Zf+eta0f) + (exp(thh)-rho)*Pf  + Ch) ...
                + exp(betf)*exp( sig*(Cf-Cf) + (1-gam)*(1-sig)*(Wf-Wf) )*( exp(EV1f)-exp(EV0f) ) - exp(Wf)*tau0f; 
% 9 eta1h
F(9)  = 1/exp(thf)*a2f*exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh)/(exp(thf)-1))^(1-exp(thf))*exp( exp(thf)*q+(exp(thf)-1)*(Zh+eta1h) + (exp(thf)-rho)*Phs + Cf) ...
                + exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) )*( exp(EV1h)-exp(EV0h) ) - exp(Wh)*tau1h;

% 10 eta0f
F(10) = 1/exp(thh)*a2h*exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf)/(exp(thh)-1))^(1-exp(thh))*exp(-exp(thh)*q+(exp(thh)-1)*(Zf+eta1f) + (exp(thh)-rho)*Pf  + Ch) ...
                + exp(betf)*exp( sig*(Cf-Cf) + (1-gam)*(1-sig)*(Wf-Wf) )*( exp(EV1f)-exp(EV0f) ) - exp(Wf)*tau1f;
% 11 Nh
F(11) = (1-normcdf(eta0h,0,sdeta))*(1-exp(Nh))+(1-normcdf(eta1h,0,sdeta))*exp(Nh) - exp(Nh);
% 12 Nf 
F(12) = (1-normcdf(eta0f,0,sdeta))*(1-exp(Nf))+(1-normcdf(eta1f,0,sdeta))*exp(Nf) - exp(Nf);
% 13 Lh
F(13) = exp(Lph) + (1-exp(Nh))*(1-normcdf(eta0h,0,sdeta))*tau0h + exp(Nh)*(1-normcdf(eta1h,0,sdeta))*tau1h - exp(Lh);
% 14 Lf
F(14) = exp(Lpf) + (1-exp(Nf))*(1-normcdf(eta0f,0,sdeta))*tau0f + exp(Nf)*(1-normcdf(eta1f,0,sdeta))*tau1f - exp(Lf);
% 15 Lph
F(15) = (th*exp(Wh)/(th-1))^(-1)*exp((1-rho)*Ph+Ch) + a2f*(exp(thf)*exp(Wh)/(exp(thf)-1))^(-1)*exp( q+(1-rho)*Phs+Cf) - exp(Lph);
% 16 Lpf
F(16) = (th*exp(Wf)/(th-1))^(-1)*exp((1-rho)*Pfs+Cf) + a2h*(exp(thh)*exp(Wf)/(exp(thh)-1))^(-1)*exp(-q+(1-rho)*Pf+Ch) - exp(Lpf);
% 17 Ph 
F(17) = (th*exp(Wh-Zh)/(th-1))^(1-th)*PsiT - exp((1-th)*Ph);
% 18 Pfs 
F(18) = (th*exp(Wf-Zf)/(th-1))^(1-th)*PsiT - exp((1-th)*Pfs);
% 19 Phs 
F(19) = exp(xif)^(1-exp(thf))*(exp(thf)*exp(Wh-Zh)/(exp(thf)-1))^(1-exp(thf))*exp(PsiXh) - exp((1-exp(thf))*(Phs+q));
% 20 Pf 
F(20) = exp(xih)^(1-exp(thh))*(exp(thh)*exp(Wf-Zf)/(exp(thh)-1))^(1-exp(thh))*exp(PsiXf) - exp((1-exp(thh))*(Pf-q));
% 21 PsiXh
F(21) = exp(Nh)*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta1h,(exp(thf)-1)*sdeta^2,sdeta)) + (1-exp(Nh))*exp(sdeta^2*(exp(thf)-1)^2/2)*(1-normcdf(eta0h,(exp(thf)-1)*sdeta^2,sdeta)) - exp(PsiXh);
% 22 PsiXf
F(22) = exp(Nf)*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta1f,(exp(thh)-1)*sdeta^2,sdeta)) + (1-exp(Nf))*exp(sdeta^2*(exp(thh)-1)^2/2)*(1-normcdf(eta0f,(exp(thh)-1)*sdeta^2,sdeta)) - exp(PsiXf);
% 23 Ch 
F(23) = exp((1-rho)*Ph) + a2h*exp((1-rho)*Pf) - 1;
% 24 Cf
F(24) = exp((1-rho)*Pfs) + a2f*exp((1-rho)*Phs) - 1;
% 25 beth
F(25) = (1-rhob)*log(bet)+fi*log(Cssh)*0 + rhob*beth-fi*Ch - sdb*ebf/2 - beth;
% 26 betf
F(26) = (1-rhob)*log(bet)+fi*log(Cssf)*0 + rhob*betf-fi*Cf + sdb*ebf/2 - betf;
% 27 V
F(27) = exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) ) - exp(V)*(1+xib*B*exp(V-YNh));
% 28 B 
F(28) = exp(Wh+Lph) + B + 1/th*exp((1-rho)*Ph+Ch) + 1/exp(thf)*a2f*exp(q+(1-rho)*Phs+Cf) - (exp(Ch) + exp(V)*B);
% 29 q
F(29) = exp(beth)*exp( sig*(Ch-Ch) + (1-gam)*(1-sig)*(Wh-Wh) )/(1+xib*B*exp(V-YNh)) ...
        -  exp(betf)*exp( q-q + sig*(Cf-Cf) + (1-gam)*(1-sig)*(Wf-Wf) )/(1-xib*B*exp(V-YNf-q));
% 30 NXY
F(30) = (exp(EXN)-exp(IMN))/exp(YNh) - NXY;
% 31 YNh
F(31) = exp( (1-rho)*Ph+Ch) + a2f*exp( q + (1-rho)*Phs + Cf ) - exp(YNh);
% 32 YNf
F(32) = exp( (1-rho)*Pfs+Cf) + a2h*exp(-q + (1-rho)*Pf + Ch ) - exp(YNf);
% 33 YRh 
F(33) = exp(YNh-Ph) - exp(YRh);
% 34 YRf 
F(34) = exp(YNf-Pfs) - exp(YRf);
% 35 EXN 
F(35) = log(a2f)+q+(1-rho)*Phs + Cf - EXN;
% 36 IMN 
F(36) = log(a2h)+(1-rho)*Pf + Ch - IMN;
% 37 Px 
F(37) = exp(Phs+q-xif+1/(exp(thf)-1)*Nh) + sdtot*etot/2 - exp(Px);
% 38 Pm
F(38) = exp(Pf-xih+1/(exp(thh)-1)*Nf) - sdtot*etot/2 - exp(Pm);
% 39 EXR
F(39) = exp(EXN-Px) - exp(EXR);
% 40 IMR
F(40) = exp(IMN-Pm) - exp(IMR);
% 41 n0h 
F(41) = 1-normcdf(eta0h,0,sdeta) - exp(n0h);
% 42 n0f
F(42) = 1-normcdf(eta0f,0,sdeta) - exp(n0f);
% 43 n1h
F(43) = normcdf(eta1h,0,sdeta) - exp(n1h);
% 44 n1f
F(44) = normcdf(eta1f,0,sdeta) - exp(n1f);
% 45 xih
F(45) = xic+0.5*xid + xia*0 - xih;
% 46 xif
F(46) = xic-0.5*xid - xif;
% 47 Zh
F(47) = rhozc*Zc + sdzc*ezc - Zc;
% 48 Zf
F(48) = rhozd*Zd + sdzd*ezd - Zd;
% 49
F(49) = Zc+0.5*Zd - Zh;
% 50
F(50) = Zc-0.5*Zd +zgap - Zf;

% 51 thh
F(51) = th*exp(zeta*q) - exp(thh);
% 52 thf
F(52) = th*exp(-zeta*q) - exp(thf);
% 53 TOT
F(53) = Pm-Px - TOT;
% 54 EXIMR
F(54) = EXR-IMR - EXIMR;
% 55 TOTa
F(55) = Pm-Px - TOTa;
% 56 tshare
F(56) = (exp(EXN)+exp(IMN))/exp(YNh) - tshare;
% 57 NXYR
F(57) = (exp(EXR)-exp(IMR))/exp(YRh) - NXYR;
% 58 WEDGEDD  
F(58) = EXIMR - (Cf - Ch) - rho*totq - WEDGEDD;
% 59 WEDGEYY  
F(59) = EXIMR - (YRf - YRh) - rho*totq - WEDGEYY;
% 60 IMPORTSHARE
F(60) = IMR - YRh - MD;
% 61 xic
F(61) = rhoxic*xic + sdci*exic/rho - xic;
% 62 xid
F(62) = rhoxid*xid + sddi*exid/rho - xid;
% 63 YSY
F(63) = YRf - YRh - YSY;
% 64 totq
F(64) = TOT+q - totq;
% 65 ipdetrend
F(65) = YRh - YRh - ipdetrend;
% 66 DSD_ADV
F(66) = Cf- Ch - DSD_ADV;
% 67 Assets to Nominal GDP
F(67) = B/exp(YNh) - By;
% 68 Nominal EXPORT IMPORT RATIO
F(68) = EXN-IMN - EXIMN;

% 69 mtotq
F(69) = totq + atotq - mtotq;

%70 mDSD_ADV
F(70) = DSD_ADV + aysy - mDSD_ADV;

%71 xia
F(71) = rhoxia*xia - xia;

%72 XMY 
F(72) = ( exp(EXR) + exp(IMR) )/exp(YRh) - exp(XMY);

% 73 mysy
F(73) = YSY+aysy - mysy;

F(74) = TOTa - mtot;

F(75) = relF - exp(q+YNf-YNh);

F=F';

