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

T = 1:1:25;
T=T';


lwidth = 1.5;
asize = 20;
fsize = 20;
pw = 8.5;
ph = 11;
lsize = 6;

%figure('Position',[20,20,750,900],'Name',...
figure('Position',[20,20,350,400],'Name',...
    '','Color','w')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pw ph]);
set(gcf, 'PaperOrientation', 'portrait')
set(0,'DefaultAxesFontSize', 40)
ZeroLine = zeros(length(T),1);


%% Figure 9 A

plot(T,oo_.irfs.EXIMR_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.EXIMR_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.EXIMR_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-3,1]);xlim([0,25]);
yticks([-3:0.5:1]);
%legend('boxoff')
%legend('z_d','\beta','\xi_d', 'Location','northwest','FontSize',20)
saveas(gcf,strcat('fig9EXIMR.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.EXIMN_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.EXIMN_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.EXIMN_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-3,1]);xlim([0,25]);
yticks([-3:0.5:1]);
%legend('boxoff')
%legend('z_d','\beta','\xi_d', 'Location','northwest','FontSize',20)
saveas(gcf,strcat('fig9EXIMN.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.q_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.q_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.q_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-1,3]);xlim([0,25]);
yticks([-1:0.5:3]);
%legend('boxoff')
%legend('z_d','\beta','\xi_d', 'Location','northwest','FontSize',20)
saveas(gcf,strcat('fig9RER.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.TOT_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.TOT_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.TOT_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-1,1]);xlim([0,25]);
yticks([-1:0.5:1]);
legend('boxoff')
legend('z_d','\beta','\xi_d', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9TOT.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.YSY_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.YSY_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.YSY_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-1,1]);xlim([0,25]);
yticks([-1:0.5:1]);
%legend('boxoff')
%legend('z_d','\beta','\xi_d', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9YSY.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.XMY_ezd(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,oo_.irfs.XMY_ebf(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.XMY_exid(1:25)*100,'-.magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-1,1]);xlim([0,25]);
yticks([-1:0.5:1]);
%legend('boxoff')
%legend('z_d','\beta','\xi_d', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9XMY.eps'),'epsc');
hold off

%% Figure 9 B
% 
plot(T,oo_.irfs.YRh_ezc(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.YRh_exic(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([0,1.5]);xlim([0,25]);
yticks([0:0.5:1.5]);
%legend('boxoff')
%legend('z_c','\xi_c', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9YRh.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.Ch_ezc(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.Ch_exic(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([0,1.5]);xlim([0,25]);
yticks([0:0.5:1.5]);
%legend('boxoff')
%legend('z_c','\xi_c', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9Ch.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.XMY_ezc(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.XMY_exic(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-0.5,2]);xlim([0,25]);
yticks([-0.5:0.5:2]);
legend('boxoff')
legend('z_c','\xi_c', 'Location','best','FontSize',20)
saveas(gcf,strcat('fig9XMYC.eps'),'epsc');
hold off
% 
plot(T,oo_.irfs.Nh_ezc(1:25)*100,'-blue','LineWidth',lwidth), hold on
plot(T,-oo_.irfs.Nh_exic(1:25)*100,'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:25])
ytickformat('%,.1f')
ylim([-0.5,2.5]);xlim([0,25]);
yticks([-0.5:0.5:2.5]);
%legend('boxoff')
%legend('z_c','\xi_c', 'Location','southeast','FontSize',20)
saveas(gcf,strcat('fig9Nh.eps'),'epsc');
hold off
