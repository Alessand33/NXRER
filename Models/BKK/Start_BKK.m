clear all;
close all;
clc;
% CHANGE LINE 30 aCeqaI  = 0 for K intensive trade and =1 common home bias
%% Using exact variance equations
%%  Parameter Values 
bet  = 0.99;             % time preference: real interest rate = 1 - 1/beta 
sig  = 2;                % CRRA coefficient 
th   = 4/5;              % Elasticity of demand: Markup = 1/theta - 1
els  = 1.5;
rho  = 1 - 1/els ;       % Elasticity of Substitution between Home and Foreign = 1/(1-rho) 

alp  = 0.36;             % capital share 
del  = 0.025;            % capital depreciation rate 

%  modified parameters 
co   = th*(1-alp);
mu   = 1/(rho-1) - 1/(th-1);
v    = (1-th)/( 1 - th*(1-alp) );
m    = (1-v)/alp;

n1   = 1 - 0;            % persistence --> persistence rate later it will be changed to stopper rate
N    = 1;                % # of exporters 
n0   = 0;                % entry rate  
sdeta  = 0.05;           % SD of firm specific shocks 

Tr   = 0.15;             % Steady state IM/Y ratio 
L    = 1/4;              % Steady state Labor supply 
a2I     = 1;  % investment imported goods share
aCeqaI  = 0; % =1 if common home bias, 0 if a2I=1 

rhoz    = 0.906;  % persistence of TFP =z 
rhozr   = 0.088; % spillover of foreign TFP =z
sdeps   = 0.007;   % SD(z shock)
coreps  = 0.0;   % correlation of z shock 
rhosdeta= 0.95;   % sd(eps) shock persistence
sdzeta = 0.00;      % sd(eps) shock volatility to make percent in IR
xi      = 1;  % iceberg cost normalized to 1
sdxi    = 0.0; % SD(xishock)
rhoxi   = 0.975; % persistence of xi shock
exi     = 1;    % shock to xi
epsh    = 1;
epsf    = 1;
zeta    = 1;   % volatility shock 

%================================================================================
%   Steady State Computation 
%================================================================================

xpar = [bet; sig; th; els; rho; alp; del; sdeta; rhoz; sdeps; coreps; co; mu; v; m; aCeqaI];  % parameters used in SteadyS 
xss0 = [Tr; L];                                                 % fixed variables used in SteadyS
save xpar;
save xss0;

%---- Initial Values indm=1

K1 = 15.3340;
a2 = 0.5;

C    = 1.06; 
Ph   = 1.56; 
Pf   = 1.59; 
gam  = 0.35;

%--- Steady State Computation 
x0= log([C;Ph;Pf;K1;a2;gam;a2I]);
if els>1.3 x0 = [0.4267;    0.5417;    0.5417;    2.8126;   -1.3343;   -1.1101;   -1.3343]; end;
if els>1 & els<1.2; x0=[6.4807; 4.4628; 4.4628; 8.5906; -1.3203;-0.9337; -1.3203]; end;
if els<1; x0 =   [-1.1897; -0.4463;    -0.4463;    0.9201;   -2.7726;    -0.9337;    -2.7726]; end

x = fsolve('SteadyS',x0,optimset('Display','iter'));
C    = exp(x(1)); 
Ph   = exp(x(2)); 
Pf   = exp(x(3)); 
K1   = exp(x(4)); 
a2   = exp(x(5)); 
gam  = exp(x(6));
a2I  = exp(x(7));

W = (1-gam)/gam*C/(1-L);
a1 = 1;
K  = K1;
Iv  = del*K;
Lp   = L;
PI  = ( Ph^(rho/(rho-1)) + a2I^(1/(1-rho))*Pf^(rho/(rho-1)) )^((rho-1)/rho);
Dh = C + PI^(1/(1-rho))*Iv;
Df = a2^(1/(1-rho))*C + a2I^(1/(1-rho))*PI^(1/(1-rho))*Iv;
L1 = co^(v/(1-th))*W^(v/(th-1))*K1^(1-v)*( ( Ph^mu*Dh + Pf^mu*Df )^v );
EX = Pf^(1/(rho-1))*Df;
YN = 1/co*W*Lp;
Y  = YN/Ph;

a2^(1/(1-rho))
a2I^(1/(1-rho))

%  Additional Variables from SteadyS Solution
NXY = 0;
Phs = Pf;
Pfs = Ph;
IM  = EX;
Q   = 1;
Z   = 1;

tw = a2^(1/(1-rho))*C/Df;
% Reparameterization
a2 = a2^(1/(1-rho)); 
a2I = a2I^(1/(1-rho)); 
n1 = 1-n1; % taking n1 as the stopper rate 


%% Dynare Program for Innovation

m_sdeta = log(sdeta);  % mean of sd(eps) in logarithm

xpar=[gam sig co v th m a1 a2 mu m_sdeta del rho rhoz sdeps coreps bet alp rhosdeta sdzeta a2I sdxi rhoxi rhozr tw];
xss =[Y Y C C Iv Iv L L Dh Dh    Df Df Ph Pfs Phs Pf PI PI W W   Q ... 
     L1 L1 K K EX IM exp(NXY) exp(NXY) Q exp(0) exp(0) exp(EX+IM)...
     K1 K1 Z Z sdeta sdeta epsh epsf zeta zeta xi xi exi exi]; 
xss = log(xss);

save xpar;
save xss;

dynare dynare_BKK;

tsize = 14;
legsize = 12;
fsize = 12;
asize = 12;

%The element oo_.autocorr{i}(k,l) is equal to the correlation between 
%$y^k_t$ and $y^l_{t-i}$, where $y^k$ (resp. $y^l$) is the $k$-th (resp. $l$-th) endogenous variable in the declaration order.
%Qt-8,nxt Qt+1,nxt
%Qt,EXDDt-k
exmddrercorr = [oo_.autocorr{12}(2,5); oo_.autocorr{11}(2,5); oo_.autocorr{10}(2,5); oo_.autocorr{9}(2,5); oo_.autocorr{8}(2,5); oo_.autocorr{7}(2,5); oo_.autocorr{6}(2,5); oo_.autocorr{5}(2,5); oo_.autocorr{4}(2,5); oo_.autocorr{3}(2,5);
    oo_.autocorr{2}(2,5);oo_.autocorr{1}(2,5);oo_.var(2,5)/(oo_.var(2,2)*oo_.var(5,5))^0.5;
    oo_.autocorr{1}(5,2); oo_.autocorr{2}(5,2); oo_.autocorr{3}(5,2); oo_.autocorr{4}(5,2); oo_.autocorr{5}(5,2); oo_.autocorr{6}(5,2); oo_.autocorr{7}(5,2); oo_.autocorr{8}(5,2); oo_.autocorr{9}(5,2); oo_.autocorr{10}(5,2); oo_.autocorr{11}(5,2); oo_.autocorr{12}(5,2)];
nxyrercorr = [oo_.autocorr{8}(2,7); oo_.autocorr{7}(2,7);
    oo_.autocorr{6}(2,7); oo_.autocorr{5}(2,7); oo_.autocorr{4}(2,7); oo_.autocorr{3}(2,7); oo_.autocorr{2}(2,7); oo_.autocorr{1}(2,7); 
    oo_.var(2,7)/(oo_.var(2,2)*oo_.var(7,7))^0.5;
    oo_.autocorr{1}(7,2); oo_.autocorr{2}(7,2); oo_.autocorr{3}(7,2); oo_.autocorr{4}(7,2); oo_.autocorr{5}(7,2); oo_.autocorr{6}(7,2); oo_.autocorr{7}(7,2); oo_.autocorr{8}(7,2)];
nxyyrercorr = [oo_.autocorr{12}(2,8); oo_.autocorr{11}(2,8); oo_.autocorr{10}(2,8); oo_.autocorr{9}(2,8); oo_.autocorr{8}(2,8); oo_.autocorr{7}(2,8); oo_.autocorr{6}(2,8); oo_.autocorr{5}(2,8); oo_.autocorr{4}(2,8); oo_.autocorr{3}(2,8); oo_.autocorr{2}(2,8);
    oo_.autocorr{1}(2,8); oo_.var(8,2)/(oo_.var(2,2)*oo_.var(8,8))^0.5; oo_.autocorr{1}(8,2); oo_.autocorr{2}(8,2); oo_.autocorr{3}(8,2); oo_.autocorr{4}(8,2);
    oo_.autocorr{5}(8,2); oo_.autocorr{6}(8,2); oo_.autocorr{7}(8,2); oo_.autocorr{8}(8,2); oo_.autocorr{9}(8,2); oo_.autocorr{10}(8,2); oo_.autocorr{11}(8,2); oo_.autocorr{12}(8,2)];
exmrercorr = [oo_.autocorr{12}(2,9); oo_.autocorr{11}(2,9); oo_.autocorr{10}(2,9); oo_.autocorr{9}(2,9); oo_.autocorr{8}(2,9); oo_.autocorr{7}(2,9); oo_.autocorr{6}(2,9); oo_.autocorr{5}(2,9); oo_.autocorr{4}(2,9); oo_.autocorr{3}(2,9); oo_.autocorr{2}(2,9);
    oo_.autocorr{1}(2,9); oo_.var(9,2)/(oo_.var(2,2)*oo_.var(9,9))^0.5; oo_.autocorr{1}(9,2); oo_.autocorr{2}(9,2); oo_.autocorr{3}(9,2); oo_.autocorr{4}(9,2);
    oo_.autocorr{5}(9,2); oo_.autocorr{6}(9,2); oo_.autocorr{7}(9,2); oo_.autocorr{8}(9,2); oo_.autocorr{9}(9,2); oo_.autocorr{10}(9,2); oo_.autocorr{11}(9,2); oo_.autocorr{12}(9,2)];
exmysdrercorr = [oo_.autocorr{12}(2,18); oo_.autocorr{11}(2,18); oo_.autocorr{10}(2,18); oo_.autocorr{9}(2,18); oo_.autocorr{8}(2,18); oo_.autocorr{7}(2,18); oo_.autocorr{6}(2,18); oo_.autocorr{5}(2,18); oo_.autocorr{4}(2,18); oo_.autocorr{3}(2,18);
    oo_.autocorr{2}(2,18);oo_.autocorr{1}(2,18);oo_.var(2,18)/(oo_.var(2,2)*oo_.var(18,18))^0.5;
    oo_.autocorr{1}(18,2); oo_.autocorr{2}(18,2); oo_.autocorr{3}(18,2); oo_.autocorr{4}(18,2); oo_.autocorr{5}(18,2); oo_.autocorr{6}(18,2); oo_.autocorr{7}(18,2); oo_.autocorr{8}(18,2); oo_.autocorr{9}(18,2); oo_.autocorr{10}(18,2); oo_.autocorr{11}(18,2); oo_.autocorr{12}(18,2)];
exmddtw_rercorr = [oo_.autocorr{12}(2,19); oo_.autocorr{11}(2,19); oo_.autocorr{10}(2,19); oo_.autocorr{9}(2,19); oo_.autocorr{8}(2,19); oo_.autocorr{7}(2,19); oo_.autocorr{6}(2,19); oo_.autocorr{5}(2,19); oo_.autocorr{4}(2,19); oo_.autocorr{3}(2,19);
    oo_.autocorr{2}(2,19);oo_.autocorr{1}(2,19);oo_.var(2,19)/(oo_.var(2,2)*oo_.var(19,19))^0.5;
    oo_.autocorr{1}(19,2); oo_.autocorr{2}(19,2); oo_.autocorr{3}(19,2); oo_.autocorr{4}(19,2); oo_.autocorr{5}(19,2); oo_.autocorr{6}(19,2); oo_.autocorr{7}(19,2); oo_.autocorr{8}(19,2); oo_.autocorr{9}(19,2); oo_.autocorr{10}(19,2); oo_.autocorr{11}(19,2); oo_.autocorr{12}(19,2)];

if els==1.5;
    if aCeqaI ==0 
    xlswrite('TheoryPredInv.xls',[exmrercorr exmddrercorr nxyyrercorr exmysdrercorr exmddtw_rercorr],'B2:F26');
    elseif aCeqaI==1 
    xlswrite('TheoryPred.xls',[exmrercorr exmddrercorr nxyyrercorr exmysdrercorr exmddtw_rercorr],'B2:F26');
    end;
end;
figure 
%plot((-12:12)',nxyrercorr, '--r', 'LineWidth',2), hold on
plot((-12:12)',exmrercorr, '-.ok', 'LineWidth',2), hold on
plot((-12:12)',exmddrercorr, '-b', 'LineWidth',2), hold on
plot((-12:12)',nxyyrercorr, '-.+k', 'LineWidth',2), hold on
plot((-12:12)',exmysdrercorr, '--r', 'LineWidth',2), hold on
plot((-12:12)',exmddtw_rercorr, '--r', 'LineWidth',5), hold on
legend('EXMRt+k','EXDsDt+k','EXYsYt+k','EXDsYt+k')
legend('EXMRt+k','EXDsDt+k','EXYsYt+k','EXDsYt+k','location','southwest', 'AutoUpdate','off')
     legend('boxoff')
title({'Lag and Lead between Net trade flows and RER';'Theory '},'FontSize',tsize);
fxlabel=xlabel('quarter');
% set(fxlabel,'FontSize',fsize)
% set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca,'XTick',-12:4:12)
ylim([-.5, 1])
xlim([-12 12])
plot([-12 12],[0 0],':k', 'LineWidth',0.25)
%print -depsc -tiff 'NXRERTheory.eps'


% % % 	xcnxrip	xcnxr	xcnxy	xchpnxy
 num = xlsread('XcorrUS.xlsx',1,'A2:E42');
figure 
plot(num(13:29,1),num(13:29,4), '--r', 'LineWidth',2), hold on
%plot(num(13:29,1),num(13:29,3), '-b', 'LineWidth',2), hold on
plot(num(13:29,1),num(13:29,2), '-.+k', 'LineWidth',2), hold on
plot(num(13:29,1),num(13:29,3), '-.ok', 'LineWidth',2), hold on
legend('NXYt+k','NXYYt+k','EXMRt+k','location','Northwest','AutoUpdate','off')
     legend('boxoff')
title({'Lag and Lead between Net trade flows and RER';'Data'},'FontSize',tsize);
fxlabel=xlabel('quarter');
% set(fxlabel,'FontSize',fsize)
% set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca,'XTick',-8:4:8)
ylim([-.5, 1])
plot([-8 8],[0 0],':k', 'LineWidth',0.25)
% print -depsc -tiff 'NXRERData.eps'
% 
% % 	xcnxrip	xcnxr	xcnxy	xchpnxy
figure 
plot(num(13:29,1),num(13:29,3), '-.r', 'LineWidth',2), hold on
plot((-8:8)',exmrercorr(5:21), '-.ok', 'LineWidth',2), hold on
legend('Data','Theory','location','Southwest')
     legend('boxoff')
title({'Lag and Lead between ln(EX/M)_t_+_k and RER_t';'Data'},'FontSize',tsize);
fxlabel=xlabel('k (quarter)');
% set(fxlabel,'FontSize',fsize)
% set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca,'XTick',-8:4:8)
ylim([-1, 1])
plot([-8 8],[0 0],':k', 'LineWidth',0.25), hold on
plot([0 0],[-1 1],':k', 'LineWidth',0.25)
% %print -depsc -tiff 'NXRERDatavsTheory.eps'
% 
figure 
plot(num(13:29,1),num(13:29,2), '-.r', 'LineWidth',2), hold on
plot((-8:8)',nxyyrercorr(5:21), '-.ok', 'LineWidth',2), hold on
legend('Data','Theory','location','Southwest')
     legend('boxoff')
title({'Lag and Lead between ln((EX/IP*)/(M/IP))_t_+_k and RER_t';'Data vs Theory'},'FontSize',tsize);
fxlabel=xlabel('k (quarter)');
% set(fxlabel,'FontSize',fsize)
% set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca,'XTick',-8:4:8)
ylim([-1, 1])
plot([-8 8],[0 0],':k', 'LineWidth',0.25), hold on
plot([0 0],[-1 1],':k', 'LineWidth',0.25)
% %print -depsc -tiff 'NXYYRERDatavsTheory.eps'



filename = 'dynare_BKK';

eval(['delete ' filename '*.m']);
eval(['delete ' filename '*.mat']);
eval(['delete ' filename '*.swp']);
eval(['delete ' filename '*.log']);
eval(['delete ' filename '*.eps']);
eval(['delete ' filename '*.pdf']);
eval(['delete ' filename '*.fig']);

disp('Steady State Statistics:')
fprintf('\n');

fprintf('Share imported goods                                     = %9.3f %9.3f \n', [0.2,          0.265]);
%report statistics

fprintf('\n');
fprintf('\n');
disp('First column: Model, Second: Data since 1995q1')
fprintf('\n');
fprintf('\n');

disp('Standard Deviations:')

stddev = sqrt(diag(oo_.var));   % order in which arranged in stoch_simul
fprintf('Based on US data from 1980q1 to 2014q3\n');
%Yh Q EX IM EXDD tot NXY NXRYY EXR tr Yf Ch Ih Cf If;
%1   2 3  4  5    6   7   8    9   10 11 12 13 14 15
fprintf('Absolute Volatility                                               Data    Theory\n');
fprintf('Production,  Y                                             = %9.2f %9.2f \n', [1.5, stddev(1)*100] );
fprintf('RER, Q                                                     = %9.2f %9.2f \n', [10.2, stddev(2)*100]);
fprintf('Volatility relative to RER\n');
fprintf('EXMR/Q, EXR                                                = %9.2f %9.2f \n', [1.34, stddev(9)/stddev(2)]);
fprintf('TOT/Q, tot                                                 = %9.2f %9.2f \n', [0.32, stddev(6)/stddev(2)]);
fprintf('Volatility relative to Y\n');
fprintf('C                                                          = %9.2f %9.2f \n', [0, stddev(12)/stddev(1)]);
fprintf('I                                                          = %9.2f %9.2f \n', [0, stddev(13)/stddev(1)]);
fprintf('Trade/Y, tr                                                = %9.2f %9.2f \n', [5.9, stddev(10)*100]);

fprintf('\n');
disp('Corelation')
fprintf('EXMR(-8),RER                                               = %9.2f %9.2f \n', [num(13,3), oo_.autocorr{8}(2,9)]);
fprintf('EXMR,RER                                                   = %9.2f %9.2f \n', [num(21,3), oo_.var(9,2)/(stddev(2)*stddev(9))]);
fprintf('EXMR(+6),RER                                               = %9.2f %9.2f \n', [num(27,3), oo_.autocorr{6}(9,2)]);
fprintf('RER,TOT                                                    = %9.2f %9.2f \n', [0.75, oo_.var(6,2)/(stddev(2)*stddev(6))]);
fprintf('Y and Ys                                                   = %9.2f %9.2f \n', [0, oo_.var(1,11)/(stddev(1)*stddev(11))]);

%fprintf('\n');
disp('AutoCorrelations:');
%fprintf('\n');
% fprintf('Production,  M                                           = %9.2f %9.2f \n', [oo_.autocorr{1}(6,6), 0.91] );
fprintf('EXMR, EX/M                                                 = %9.2f %9.2f \n', [0.9831, oo_.autocorr{1}(9,9)]);
fprintf('RER                                                        = %9.2f %9.2f \n', [0.9721, oo_.autocorr{1}(2,2)]);
fprintf('TOT                                                        = %9.2f %9.2f \n', [0.91, oo_.autocorr{1}(6,6)]);



%fprintf('Abs.Vol\n');
disp('  ');
fprintf('%9.2f\n', [stddev(2)*100]);
%fprintf('Vol.rel.Q\n');
disp('  ');
fprintf('%9.2f\n', [stddev(9)/stddev(2)]);
fprintf('%9.2f\n', [stddev(6)/stddev(2)]);
%fprintf('Vol.rel.Y\n');
disp('  ');
fprintf('%9.2f\n', [stddev(12)/stddev(1)]);
fprintf('%9.2f\n', [stddev(13)/stddev(1)]);
disp('  ');
fprintf('%9.2f\n', [oo_.autocorr{8}(2,9)]);
fprintf('%9.2f\n', [oo_.var(9,2)/(stddev(2)*stddev(9))]);
fprintf('%9.2f\n', [oo_.autocorr{6}(9,2)]);
fprintf('%9.2f\n', [oo_.var(6,2)/(stddev(2)*stddev(6))]);
fprintf('%9.2f\n', [oo_.var(1,11)/(stddev(1)*stddev(11))]);

%fprintf('\n');
%disp('AutoCorrelations:');
disp('  ');
%fprintf('\n');
fprintf('%9.2f\n',[oo_.autocorr{1}(9,9)]);
fprintf('%9.2f\n',[oo_.autocorr{1}(2,2)]);
fprintf('%9.2f\n', [oo_.autocorr{1}(6,6)]);

