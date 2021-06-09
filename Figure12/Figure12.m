%% Figure 12
clear all;
close all;
clc;

Data0= xlsread('Sorted.xlsx',2);
sizedata=size(Data0);
T=1980:0.25:2014.5;
T=T';

lwidth = 3;
asize = 20;
fsize = 20;
pw = 8.5;
ph = 11;
lsize = 6;

%figure('Position',[20,20,750,900],'Name',...
figure('Position',[20,20,550,425],'Name',...
    '','Color','w')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pw ph]);
set(gcf, 'PaperOrientation', 'portrait')
set(0,'DefaultAxesFontSize', 40)
ZeroLine = zeros(length(T),1);

% RER
plot(T,Data0(:,1),'-blue','LineWidth',lwidth), hold on
plot(T,Data0(:,2),'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',1)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-30,30]);xlim([1980,2015]);
yticks([-30:10:30]);
ylabel('Percent') 
saveas(gcf,strcat('fig12RER.eps'),'epsc');
hold off

% TOT
plot(T,Data0(:,3),'-blue','LineWidth',lwidth), hold on
plot(T,Data0(:,4),'--red','LineWidth',lwidth), hold on
plot(T,Data0(:,5),':black','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',1)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-40,40]);xlim([1980,2015]);
yticks([-40:20:40]);
ylabel('Percent') 
legend('boxoff')
legend('Data','Model','Measurement', 'Location','northwest','FontSize',16)
saveas(gcf,strcat('fig12TOT.eps'),'epsc');
hold off

% EXIMN
plot(T,Data0(:,6),'-blue','LineWidth',lwidth), hold on
plot(T,Data0(:,7),'--red','LineWidth',lwidth), hold on
plot(T,ZeroLine,':black','LineWidth',1)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-60,20]);xlim([1980,2015]);
yticks([-60:20:20]);
ylabel('Percent') 
saveas(gcf,strcat('fig12EXIMN.eps'),'epsc');
hold off