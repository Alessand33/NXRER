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


%%  
load dynare_Benchmark_results;

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


% % Loading shocks
ex_=[];
for shock_iter=1:6
    ex_=[ex_ oo_.SmoothedShocks.(deblank(x.M_.exo_names(shock_iter,:)))];
end
% save ex_;
% % use shocks only starting at t=2 due to t=1 being initial condition
ex1_ = ex_(1:end,:);
% 
% finding initial values of variables
y0=[];
for endo_iter=1:74
	y0 = [y0;
         oo_.SmoothedVariables.(deblank(M_.endo_names(endo_iter,:)))(1)];
end;
% save y0;
% 
ysteady_state = oo_.steady_state;
% save ysteady_state;

% initial_condition_states = repmat(oo_.dr.ys,1,x.M_.maximum_lag);

% Bug correction 
y0(13)= y0(13)+log(Nh);
y0(14)= y0(14)+log(Nf);
y0(15)= y0(15)+log(beth);
y0(16)= y0(16)+log(betf);

 ex_all = ex_(2:end,:);
 dr = oo_.dr;
 % y0 = initial_condition_states;
 % 1. All Shocks + Initial
 iorder = 1;
 ex_tmp = ex_all*0;
 yd_tmp = y0; 
 
 % yd_tmp = ysteady_state;
 Actual=simult_(yd_tmp,dr,ex_tmp,iorder);
 Actual=Actual';
  
 
S_xid = [Actual(:,60) oo_.SmoothedVariables.xid];
S_xic = [Actual(:,59) oo_.SmoothedVariables.xic];
S_zd =[Actual(:,66) oo_.SmoothedVariables.Zd]; 
S_zc =[Actual(:,65) oo_.SmoothedVariables.Zc]; 
S_Nh = exp([Actual(:,13) oo_.SmoothedVariables.Nh+log(Nh)]);  % Bug correction
S_Nf = exp([Actual(:,14) oo_.SmoothedVariables.Nf+log(Nf)]);  % Bug correction
S_beth = [Actual(:,15)-log(beth) oo_.SmoothedVariables.beth]; % Bug correction
S_betf = [Actual(:,16)-log(betf) oo_.SmoothedVariables.betf]; % Bug correction
S_betd = S_beth-S_betf;
T = 1980:0.25:2014.5;
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


%% Figure 8
NhData= xlsread('NhData.xlsx',1);

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_xid(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_xid(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-1,0.2]);xlim([1980,2015]);
yticks([-1:0.2:0.2]);
saveas(gcf,strcat('fig8xid.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_xic(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_xic(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.1,0.3]);xlim([1980,2015]);
yticks([-0.1:0.1:0.3]);
saveas(gcf,strcat('fig8xic.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_beth(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_beth(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-0.0065,0.005]);xlim([1980,2015]);
yticks([-0.006:0.002:0.004]);
saveas(gcf,strcat('fig8beth.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_zd(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_zd(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.2f')
ylim([-0.5,0.0]);xlim([1980,2015]);
yticks([-0.5:0.1:0.0]);
saveas(gcf,strcat('fig8zd.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_zc(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_zc(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.03,0.3]);xlim([1980,2015]);
yticks([-0.0:0.1:0.3]);
saveas(gcf,strcat('fig8zc.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_betd(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_betd(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-0.005,0.02]);xlim([1980,2015]);
yticks([-0.005:0.005:0.02]);
saveas(gcf,strcat('fig8betd.eps'),'epsc');
hold off

plot(T,S_Nh(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_Nh(:,2),'-blue','LineWidth',lwidth), hold on
plot(NhData(:,1),NhData(:,2),'-o magenta','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',0.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
legend('boxoff')
legend('Initial','Model','Data', 'Location','northwest','FontSize',20)
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([0.00,0.4]);xlim([1980,2015]);
yticks([0.0:0.1:0.4]);
saveas(gcf,strcat('fig8Nh.eps'),'epsc');
hold off


plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,S_Nf(:,1),'--red','LineWidth',lwidth), hold on
plot(T,S_Nf(:,2),'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.2f')
ylim([0.00,0.25]);xlim([1980,2015]);
yticks([0.0:0.05:0.25]);
saveas(gcf,strcat('fig8Nf.eps'),'epsc');
hold off


