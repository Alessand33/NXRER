clear all;
close all;
clc;

CorrData= xlsread('Figure3.xlsx',3);

lwidth = 3;
asize = 15;
fsize = 15;
pw = 8.5;
ph = 11;
lsize = 6;

%figure('Position',[20,20,750,900],'Name',...
figure('Position',[20,20,550,425],'Name',...
    '','Color','w')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pw ph]);
set(gcf, 'PaperOrientation', 'portrait')
set(0,'DefaultAxesFontSize', 12)
ZeroLine = zeros(length(CorrData),1);

plot(CorrData(:,1),CorrData(:,2),'-blue','LineWidth',lwidth), hold on
plot(CorrData(:,1),CorrData(:,3),'--red','LineWidth',lwidth), hold on
plot(CorrData(:,1),CorrData(:,4),':black','LineWidth',lwidth)
set(gca, 'YGrid', 'on', 'XGrid', 'on')
set(gca,'FontSize',asize)
set(gca, 'Box', 'On');
xticks([-12:4:12])
ytickformat('%,.1f')
ylim([-0.6,1.1]);xlim([-12,12]);
yticks([-0.5:0.5:1]);
legend('boxoff')
legend('Data','BKK','BKK Invest. Intensive Trades', 'Location','northwest','FontSize',15)
xlabel('Quarter(k)') 
ylabel('Correlation') 
saveas(gcf,strcat('fig3A.eps'),'epsc');
hold off

plot(CorrData(:,1),CorrData(:,5),'-blue','LineWidth',lwidth), hold on
plot(CorrData(:,1),CorrData(:,6),'--red','LineWidth',lwidth), hold on
plot(CorrData(:,1),CorrData(:,7),':black','LineWidth',lwidth)
set(gca, 'YGrid', 'on', 'XGrid', 'on')
set(gca,'FontSize',asize)
set(gca, 'Box', 'On');
xticks([-12:4:12])
ytickformat('%,.1f')
ylim([-0.6,1.1]);xlim([-12,12]);
yticks([-0.5:0.5:1]);
legend('boxoff')
legend('Data','BKK','BKK Invest. Intensive Trades', 'Location','southeast','FontSize',15)
xlabel('Quarter(k)') 
ylabel('Correlation') 
saveas(gcf,strcat('fig3B.eps'),'epsc');
hold off
