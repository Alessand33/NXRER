
clear all;
close all;
clc;
FigData=[];

x = load('dynare_Asymmetric00_IR_results.mat');
FigData = x.oo_.irfs.EXIMN_exid'*100;

x = load('dynare_Asymmetric20_IR_results.mat');
FigData(:,2) = x.oo_.irfs.EXIMN_exid'*100;

x = load('dynare_Asymmetric40_IR_results.mat');
FigData(:,3) = x.oo_.irfs.EXIMN_exid'*100;
FigData(:,4) = x.oo_.irfs.EXIMN_exic'*100;

T=1:1:40;
T=T';

lwidth = 3;
asize = 11;
fsize = 11;
pw = 8.5;
ph = 11;
lsize = 6;
ZeroLine = zeros(length(T),1);

%figure('Position',[20,20,750,900],'Name',...
figure('Position',[20,20,550*1,425*1],'Name',...
    '','Color','w')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pw ph]);
set(gcf, 'PaperOrientation', 'portrait')
set(0,'DefaultAxesFontSize', 40)
ZeroLine = zeros(length(T),1);


% EXIMN
plot(T,FigData(1:40,1),'-black','LineWidth',lwidth), hold on
plot(T,FigData(1:40,2),'--blue','LineWidth',lwidth), hold on
plot(T,FigData(1:40,3),':red','LineWidth',lwidth), hold on
plot(T,FigData(1:40,4),'--magenta','LineWidth',lwidth, 'marker','o'), hold on
plot(T,ZeroLine,':black','LineWidth',1.5)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([0:10:40])
ytickformat('%,.1f')
ylim([-0.2,1]);xlim([0,40]);
yticks([-0.2:0.2:1]);
ylabel('Percent') 
legend('boxoff')
legend('Symmetric \xi_d shock','Asymmetric 20% \xi_d shock','Asymmetric 40% \xi_d shock','Asymmetric 40% \xi_c shock', 'Location','northeast','FontSize',11)
saveas(gcf,strcat('figA3.eps'),'epsc');
hold off
