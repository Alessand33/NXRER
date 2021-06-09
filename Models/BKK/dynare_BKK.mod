%%-----------------------------------------------------------------------------------------------------------------------------------------------
var  Yh, Yf, Ch, Cf, Ivh, Ivf, Lh, Lf, Dh, Dfs, Dhs, Df, Ph, Pfs, Phs, Pf, PI, PIs, Wh, Wf, Q
     L1h, L1f, Kh, Kf, EX, IM, NXY, EXR, tot, EXDD, NXRYY, tr
     K1h, K1f, Zh, Zf, sdetah, sdetaf, xih, xif, NXRYsD, D, Ds, DTW, DTWs, EXDD_TW;
varexo epsh, epsf, zetah, zetaf, exih, exif;

load xss;

%%-----------------------------------------------------------------------------------------------------------------------------------------------
parameters gam, sig, co, v, th, m, a1, a2, mu, m_sdeta, del, rho, rhoz, sdeps, coreps, bet, alp, rhosdeta, sdzeta, a2I, sdxi, rhoxi, rhozr, tw;

load xpar;
gam   = xpar(1);
sig   = xpar(2);
co    = xpar(3);
v     = xpar(4);
th    = xpar(5);
m     = xpar(6);
a1    = xpar(7);
a2    = xpar(8);
mu    = xpar(9);
m_sdeta = xpar(10);
del   = xpar(11);
rho   = xpar(12);
rhoz  = xpar(13);
sdeps = xpar(14);
coreps = xpar(15);
bet   = xpar(16);
alp   = xpar(17);
rhosdeta = xpar(18);
sdzeta = xpar(19);
a2I  = xpar(20);
sdxi  = xpar(21);
rhoxi  = xpar(22);
rhozr = xpar(23);
tw = xpar(24);

%%-----------------------------------------------------------------------------------------------------------------------------------------------
model; 
 
(1-gam)/gam*exp(Ch-Wh) = ( 1 - exp(Lh) );

(1-gam)/gam*exp(Cf-Wf) = ( 1 - exp(Lf) ); 

sig*( Ch-Cf ) + (1-gam)*(1-sig)*( Wh-Wf ) = Q;

co^(v/(1-th))*exp( v/(th-1)*Wh + m*Zh + (1-v)*K1h(-1) )*(  ( exp(mu*Ph+Dh) + exp( 1/(1-th)*Q + mu*Phs+Dhs ) )^v ) = exp(L1h);

co^(v/(1-th))*exp( v/(th-1)*Wf + m*Zf + (1-v)*K1f(-1) )*(  ( exp(mu*Pfs+Dfs) + exp( -1/(1-th)*Q + mu*Pf+Df ) )^v ) = exp(L1f);

exp(Ch) + exp(1/(1-rho)*PI+Ivh) = exp(Dh);

exp(Cf) + exp(1/(1-rho)*PIs+Ivf) = exp(Dfs);

a2*exp(rho/(rho-1)*xih)*exp(Ch) + a2I*exp(rho/(rho-1)*xih)*exp(1/(1-rho)*PI+Ivh) = exp(Df);

a2*exp(rho/(rho-1)*xif)*exp(Cf) + a2I*exp(rho/(rho-1)*xif)*exp(1/(1-rho)*PIs+Ivf) = exp(Dhs);

exp(L1h) = exp(Lh);

exp(L1f) = exp(Lf);

exp(K1h) = exp(Kh);

exp(K1f) = exp(Kf);

exp(Kh) - (1-del)*exp(Kh(-1)) = exp(Ivh);

exp(Kf) - (1-del)*exp(Kf(-1)) = exp(Ivf);

exp( th/(th-1)*(Q+Phs) ) = exp(th/(th-1)*Ph);

exp( th/(th-1)*(-Q+Pf) ) = exp(th/(th-1)*Pfs);

co^(v/(1-th)-1)*exp( th/(1-th)*Q + (1+v/(th-1))*Wh + m*Zh )*( exp(mu*Ph+Dh) + exp( 1/(1-th)*Q + mu*Phs + Dhs ) )^(v-1)
   *( exp((1-v)*K1h(-1)) ) = exp( th/(th-1)*Phs);
            
co^(v/(1-th)-1)*exp(- th/(1-th)*Q + (1+v/(th-1))*Wf + m*Zf )*( exp(mu*Pfs+Dfs) + exp(- 1/(1-th)*Q + mu*Pf + Df ) )^(v-1)
   *( exp((1-v)*K1f(-1)) ) = exp( th/(th-1)*Pf);
            
exp( rho/(rho-1)*Ph ) + a2*exp(rho/(rho-1)*xih)*exp( rho/(rho-1)*Pf ) = 1;

exp( rho/(rho-1)*Pfs ) + a2*exp(rho/(rho-1)*xif)*exp( rho/(rho-1)*Phs ) = 1;

exp( rho/(rho-1)*Ph ) + a2I*exp(rho/(rho-1)*xih)*exp( rho/(rho-1)*Pf ) = exp( rho/(rho-1)*PI );

exp( rho/(rho-1)*Pfs ) + a2I*exp(rho/(rho-1)*xif)*exp( rho/(rho-1)*Phs ) = exp( rho/(rho-1)*PIs );

bet*exp( sig*(Ch-Ch(+1)) + (1-gam)*(1-sig)*(Wh-Wh(+1)) )*( alp/(1-alp)*exp( Wh(+1)+L1h(+1)-PI(+1)-K1h ) +  1-del ) = 1;

bet*exp( sig*(Cf-Cf(+1)) + (1-gam)*(1-sig)*(Wf-Wf(+1)) )*( alp/(1-alp)*exp( Wf(+1)+L1f(+1)-PIs(+1)-K1f ) +  1-del ) = 1;

Zh = rhoz*Zh(-1) +rhozr*Zf(-1)+ sdeps*epsh;

Zf = rhoz*Zf(-1) +rhozr*Zh(-1)+ sdeps*epsf;

exp( 1/(rho-1)*Phs + Dhs ) = exp(EX);  

exp( 1/(rho-1)*Pf + Df ) = exp(IM);   

( exp(Q+Phs+EX) - exp(Pf+IM) )/( exp(Ph+Yh) ) = NXY; 

exp(Wh-Ph)*( exp(L1h) ) = (th*(1-alp))*exp(Yh);

exp(Wf-Pfs)*( exp(L1f) ) = (th*(1-alp))*exp(Yf);

sdetah = (1-rhosdeta)*m_sdeta + rhosdeta*sdetah(-1) + sdzeta*zetah;

sdetaf = (1-rhosdeta)*m_sdeta + rhosdeta*sdetaf(-1) + sdzeta*zetaf;

xih = rhoxi*xih(-1) + sdxi*exih;

xif = rhoxi*xif(-1) + sdxi*exif;

EXR = EX - IM;

tot = Pf - Phs;

EXDD = EXR - (Ds - D);

NXRYY = EXR - (Yf - Yh);

exp(tr) = exp(EX)+exp(IM);

NXRYsD = EXR - (Yf - D);

exp(D) = exp(Ch) + exp(Ivh);

exp(Ds) = exp(Cf) + exp(Ivf);

DTW = tw*Ch + (1-tw)*Ivh;

DTWs = tw*Cf + (1-tw)*Ivf;

EXDD_TW = EXR - (DTWs - DTW);

end;
%---------------------------------------------
%  -------------  END OF MODEL   -------------
%---------------------------------------------


initval;

Yh    = xss(1);
Yf    = xss(2);
Ch    = xss(3);
Cf    = xss(4);
Ivh   = xss(5);
Ivf   = xss(6);
Lh    = xss(7);
Lf    = xss(8);
Dh    = xss(9);
Dfs   = xss(10);
Dhs   = xss(11);
Df    = xss(12);
Ph    = xss(13);
Pfs   = xss(14);
Phs   = xss(15);
Pf    = xss(16);
PI    = xss(17);
PIs   = xss(18);
Wh    = xss(19);
Wf    = xss(20);
Q     = xss(21);
L1h   = xss(22);
L1f   = xss(23);
Kh    = xss(24);
Kf    = xss(25);
EX    = xss(26);
IM    = xss(27);
NXY   = xss(28);
EXR   = xss(29);
tot   = xss(30);
EXDD  = xss(31);
NXRYY = xss(32);
tr    =xss(33);
K1h   = xss(34);
K1f   = xss(35);
Zh    = xss(36);
Zf    = xss(37);
sdetah =  xss(38);
sdetaf =  xss(39);
epsh  = xss(40);
epsf  = xss(41);
zetah  =  xss(42);
zetaf  =  xss(43);
xih  =  xss(44);
xif  =  xss(45);
exih =  xss(46);
exif =  xss(47);
NXRYsD = EXR - (Yf - D);

D = ln(exp(Ch) + exp(Ivh));

Ds = ln(exp(Cf) + exp(Ivf));

DTW = tw*Ch + (1-tw)*Ivh;

DTWs = DTW;

EXDD_TW = EXR - (DTWs - DTW);

end;

%steady;

%-----------------------
%----    SHOCKS   ------
%-----------------------
shocks;

var epsh = 1;  % aggregate productivity 
var epsf = 1;
var zetah = 1;  % firm level productivity 
var zetaf = 1;
var zetah,zetaf = 0;
var exih = 1;  % trade cost 
var exif = 1;
var epsh,epsf = 0;
var exih,exif = 0.0;
var exih,epsh = 0;
var exif,epsf = 0;
end;

%------------------------------
%---- GETTING INNOVATIONS  ----
%------------------------------
%stoch_simul(order=1);
%close all;

%stoch_simul(order=1,periods=1000,drop=200,nograph);
%sto_mean = oo_.mean;
%stoch_simul(order=1,hp_filter=1600,periods=10,drop=1,nograph) Yh Yf;
%stoch_simul(order=1,hp_filter=1600,periods=10,drop=1,nograph);
%stoch_simul(order=1,nograph,hp_filter=1600,ar=12) Yh Q EX IM EXDD tot NXY NXRYY EXR tr;
stoch_simul(order=1,nograph,ar=12, noprint) Yh Q EX IM EXDD tot NXY NXRYY EXR tr Yf Ch Ivh Cf Ivf Lh Lf NXRYsD EXDD_TW;
%                                            1 2  3  4  5    6   7    8    9   10 11 12 13 14  15 16 17 18     19

%===================================================================%
%%%%           DEFINE STDEV OF STRUCTURAL INNOVATIONS            %%%%
%===================================================================%

%shocks;
%var e1;
%stderr 1;
%var e2;
%stderr 1;
%var e4;
%stderr 1;
%var e5;
%stderr 1;
%end;

%===================================================================%
%%%%                     SOLVE THE MODEL                         %%%%
%===================================================================%

%steady;

%estimated_params;
%rho_1, 0.9973,normal_pdf, 0.5, 2;
%rho_2, 0.9438,normal_pdf, 0.5, 2;
%rho_4, 0.95,normal_pdf,0.5,2;
%rho_5, 0.9995,normal_pdf, 0.5, 2;
%sd_e1, 0.0021,inv_gamma_pdf, 0.01, inf; 
%sd_e2, 0.0139,inv_gamma_pdf, 0.01, inf;
%sd_e4, 0.01,inv_gamma_pdf,0.01,inf;
%sd_e5, 0.0026,inv_gamma_pdf, 0.01, inf;
%alp2, 0.96,inv_gamma_pdf, 1, inf;
%end;

//varobs Y_obs L_obs ex_obs L2_obs en_obs;
%varobs Y_obs L_obs ex_obs en_obs;

%estimation(datafile=datamod9,nobs=35,mh_replic=10000,mh_nblocks=5,mh_jscale=0.8,prefilter = 1);
