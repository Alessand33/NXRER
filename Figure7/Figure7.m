clear all;
close all;
clc;

x=load('dynare_Benchmark_results.mat', 'oo_');

EXIMR_S=x.oo_.SmoothedVariables.EXIMR; % EXIMR
XMY_S=x.oo_.SmoothedVariables.XMY; % XMY
ipdetrend_S=x.oo_.SmoothedVariables.ipdetrend; % ipdetrend
DSD_ADV_S=x.oo_.SmoothedVariables.DSD_ADV; % DSD_ADV
DSD_ADV_S= -(DSD_ADV_S-mean(DSD_ADV_S)); 
totq_S=x.oo_.SmoothedVariables.totq; % totq
totq_S=totq_S - mean(totq_S); % totq
TOT_S = x.oo_.SmoothedVariables.TOT; % totq
TOT_S=TOT_S - mean(TOT_S); % totq

ezc_S=x.oo_.SmoothedShocks.ezc; % ezc
ezd_S=x.oo_.SmoothedShocks.ezd; % ezd
ebf_S=x.oo_.SmoothedShocks.ebf; % ebf
exic_S=x.oo_.SmoothedShocks.exic; % exic
exid_S= x.oo_.SmoothedShocks.exid; % exid
etot_S=x.oo_.SmoothedShocks.etot; % etot

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


%% Figure 7.A

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,EXIMR_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.5,0.1]);xlim([1980,2015]);
saveas(gcf,strcat('fig7AEXIMR.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,XMY_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.7,0.6]);xlim([1980,2015]);
saveas(gcf,strcat('fig7AXMY.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,ipdetrend_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
yticks([-0.2:0.1:0.2]);
ylim([-0.2,0.2]);xlim([1980,2015]);
saveas(gcf,strcat('fig7AUSIP.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,DSD_ADV_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.2,0.2]);xlim([1980,2015]);
yticks([-0.2:0.1:0.2]);
saveas(gcf,strcat('fig7ADSD_ADV.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,totq_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.1f')
ylim([-0.4,0.2]);xlim([1980,2015]);
yticks([-0.4:0.1:0.2]);
saveas(gcf,strcat('fig7Atotq.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,TOT_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
yticks([-0.1:0.05:0.1]);
ytickformat('%,.2f')
ylim([-0.1,0.1]);xlim([1980,2015]);
saveas(gcf,strcat('fig7ATOT.eps'),'epsc');
hold off





%% Figure 7.B

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,ezc_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-6,2]);xlim([1980,2015]);
yticks([-6:2:2]);
saveas(gcf,strcat('fig7Bezc.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,ezd_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-6,4]);xlim([1980,2015]);
yticks([-6:2:4]);
saveas(gcf,strcat('fig7Bezd.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,ebf_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-3,4]);xlim([1980,2015]);
yticks([-3:1:4]);
saveas(gcf,strcat('fig7Bebf.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,exic_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-3,2]);xlim([1980,2015]);
yticks([-3:1:2]);
saveas(gcf,strcat('fig7Bezic.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,exid_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-2,4]);xlim([1980,2015]);
yticks([-2:2:4]);
saveas(gcf,strcat('fig7Bezid.eps'),'epsc');
hold off

plot(T,ZeroLine,':black','LineWidth',0.5), hold on
plot(T,etot_S,'-blue','LineWidth',lwidth)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
yticks([-4:2:4]);
ytickformat('%,.0f')
ylim([-4,4]);xlim([1980,2015]);
saveas(gcf,strcat('fig7Betot.eps'),'epsc');
hold off

