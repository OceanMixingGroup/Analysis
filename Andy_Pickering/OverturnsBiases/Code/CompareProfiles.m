%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareProfiles.m
%
%
% * Add resampled profiles
%
% 13 Jan 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all


cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=4;

% Load full T-chain data
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))

% Load resampled data
whcase=1
testnum=1
clear fname x_resamp
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
load(fullfile('Data/',fname))

%addpath /Users/Andy/Cruises_Research/CTDstruct/

%for whp=13900:14200%2400:2600%2000:2200%
for whp=1300%:5:2000
    
    tnow=xx2.yday(whp)
    
    xl=tnow+0.25*[-1 1]
    
    id=isin(xx2.yday,xl);
    tm=nanmean(xx2.T,2);
    tm=tm(~isnan(tm));
    
    figure(1);clf
    %ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.04, 3,1);
    
    %axes(ax(1))
    ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
    colorbar
    caxis([1 7])
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:10:end),'k')
    hold on
    plot(x_resamp.tsampall,x_resamp.zall,'w')
    
    freqline(tnow)
    colormap(jet)
    ylabel('Depth (m) ','fontsize',16)
    xlabel('Yearday  ','fontsize',16)
    %
    
    figure(2);clf
    plot(xx2.T(:,whp),xx2.z,'o-')
    hold on
    %plot(x_resamp.T(
    axis ij
    ylim([1200 2000])
    
    pause(1)
    
end

%% Slightly different version: plot only profiles corresponding to each resampled profile

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=4;

% Load full T-chain data
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))

% Load resampled data
whcase=1
testnum=1
clear fname x_resamp
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
load(fullfile('Data/',fname))

for whp=1:length(x_resamp.time)
    tnow=x_resamp.time(whp)
    xl=tnow+0.25*[-1 1]
    id=isin(xx2.yday,xl);
    
    [val,I]=nanmin(abs(xx2.yday-tnow));
    
    tm=nanmean(xx2.T,2);
    tm=tm(~isnan(tm));
    
    figure(1);clf
    %ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.04, 3,1);
    
    %axes(ax(1))
    ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
    colorbar
    caxis([1 7])
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:10:end),'k')
    hold on
    plot(x_resamp.tsampall,x_resamp.zall,'w')
    
    freqline(tnow)
    colormap(jet)
    ylabel('Depth (m) ','fontsize',16)
    xlabel('Yearday  ','fontsize',16)
    
    % Plot temperature profiles
    figure(2);clf
    subplot(121)
    plot(x_resamp.T(:,whp),x_resamp.z)
    hold on
    plot(xx2.T(:,I),xx2.z)
    axis ij
    grid on
    legend('True','Resamp')

    subplot(122)
    plot(x_resamp.Lttot(:,whp),x_resamp.z,'o')
    hold on
    plot(xx2.Lttot(:,I),xx2.z,'o')
    axis ij
    xlim([0 100])
    ylim([600 2200])
    grid on
    
    pause(1)
end
%%

