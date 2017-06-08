%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotResamp_CTDscenario.m
%
% Plot effects of resampling on T-chain, specifically for CTD scenario
% using only downcasts.
%
% See also PlotResampResults12Nov.m
%
% 18 Jan 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=1
whmoor=3;xl=[1e-9 10^(-5.5)]
%whmoor=4;xl=[1e-9 10^(-6.7)]
minOT=50
testnum=5 % 1m/s
%testnum=4 % 1m/s

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

figure(1);clf
agutwocolumn(0.8)
wysiwyg
set(gcf,'defaultaxesfontsize',15)
ezall=[];

for whcase=1:100
%    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
    load(fullfile('Data/',fname))
    %    semilogx(nanmean(x_resamp.eps,2),x_resamp.z,'.','color',0.5*[1 1 1])
    semilogx(nanmean(x_resamp.eps(:,1:2:end),2),x_resamp.z,'.','color',0.75*[1 1 1])
    hold on
    %    ezall=[ezall nanmean(x_resamp.eps,2)];
    ezall=[ezall nanmean(x_resamp.eps(:,1:2:end),2)];
end

%~~  Load T-chain data to resample
%load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

% find data within time of resampling period
idtr=isin(xx2.yday,[x_resamp.time(1) x_resamp.time(end)]);

% plot actual time-mean profile
semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)

% plot resampled time-mean profile
semilogx(nanmean(ezall,2),x_resamp.z,'linewidth',2)

% plot +/- std
semilogx(nanmean(ezall,2)-std(ezall,0,2),x_resamp.z,'--','color',0.3*[1 1 1],'linewidth',2)
semilogx(nanmean(ezall,2)+std(ezall,0,2),x_resamp.z,'--','color',0.3*[1 1 1],'linewidth',2)

axis ij
xlim(xl)
ylim([nanmin(x_resamp.z)  nanmax(x_resamp.z)])
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
xlabel('<\epsilon> / Wkg^{-1}','fontsize',16)
ylabel('Depth / m ','fontsize',16)
title(['Tchain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) ' ms^{-1}, downcast only'])
grid on
shg

if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth_CTDscenario' ])
    export_fig(fname,'-pdf')
end
%%