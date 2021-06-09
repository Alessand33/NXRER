
clear all;
close all;
clc;

x=load('dynare_Benchmark_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)
% ezc	ezd	exic	exid	ebf	etot	initial	data
Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 

Bench_exid = [Decomp_EXIMN(:,4)+Decomp_EXIMN(:,3) Decomp_EXIMR(:,4)+Decomp_EXIMR(:,3) ];

x=load('dynare_HighTrade_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)

Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 


HighTrade_exid = [Decomp_EXIMN(:,4)+Decomp_EXIMN(:,3) Decomp_EXIMR(:,4)+Decomp_EXIMR(:,3) ];

x=load('dynare_LowTrade_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)

Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 

LowTrade_exid = [Decomp_EXIMN(:,4)+Decomp_EXIMN(:,3) Decomp_EXIMR(:,4)+Decomp_EXIMR(:,3) ];


x=load('dynare_NoPTM_results.mat');

%xxx=permute(oo_.shock_decomposition(11,:,:),[3 2 1])
%(var, shock, year)

Decomp_EXIMR=permute(x.oo_.shock_decomposition(11,:,:),[3 2 1])*100  % EXIMR 
Decomp_EXIMN=permute(x.oo_.shock_decomposition(68,:,:),[3 2 1])*100  % EXIMN 

Armington_exid = [Decomp_EXIMN(:,4)+Decomp_EXIMN(:,3) Decomp_EXIMR(:,4)+Decomp_EXIMR(:,3) ];


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


% EXIMN
plot(T,Bench_exid(:,1),'-black','LineWidth',lwidth), hold on
plot(T,LowTrade_exid(:,1),'--blue','LineWidth',lwidth), hold on
plot(T,HighTrade_exid(:,1),':red','LineWidth',lwidth), hold on
plot(T,Armington_exid(:,1),'--black','marker','x','LineWidth',1.5), hold on
plot(T,ZeroLine,':black','LineWidth',1)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-15,15]);xlim([1980,2015]);
yticks([-15:5:15]);
ylabel('Percent') 
legend('boxoff')
legend('Benchmark','LowTrade','HighTrade', 'Location','northeast','FontSize',16)
saveas(gcf,strcat('figA2EXIMN.eps'),'epsc');
hold off

% EXIMR
plot(T,Bench_exid(:,2),'-black','LineWidth',lwidth), hold on
plot(T,LowTrade_exid(:,2),'--blue','LineWidth',lwidth), hold on
plot(T,HighTrade_exid(:,2),':red','LineWidth',lwidth), hold on
plot(T,Armington_exid(:,2),'--black','marker','x','LineWidth',1.5), hold on
plot(T,ZeroLine,':black','LineWidth',1)
%set(fylabel,'FontSize',fsize)
set(gca,'FontSize',asize)
set(gca, 'Box', 'Off');
xticks([1980:10:2010])
ytickformat('%,.0f')
ylim([-15,15]);xlim([1980,2015]);
yticks([-15:5:15]);
ylabel('Percent') 
%legend('boxoff')
%legend('Benchmark','Static-NoPTM','Static-PTM', 'Location','northwest','FontSize',16)
saveas(gcf,strcat('figA2EXIMR.eps'),'epsc');
hold off
