%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotBiasFactorTchain.m
%
% Plot the ratio of the resampled ensemble epsilon to the true value to see
% by what factor the resampled data is biased high.
%
% Plot for all different speeds on same figure.
%
%
% 3 Feb. 2015
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot mean depth-profiles of epsilon

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3;
minOT=50

figure(1);clf
agutwocolumn(0.6)
wysiwyg

for testnum=[3 1 7 6 2 4 5]
    
    ezall=[];
    
    for whcase=1:100
        fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
        load(fullfile('Data/',fname))
        hold on
        ezall=[ezall nanmean(x_resamp.eps,2)];
    end
    
    %~~  Load T-chain data to resample
    load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))
    % find data within time of resampling period
    idtr=isin(xx2.yday,[x_resamp.time(1) x_resamp.time(end)]);
    
    plot(nanmean(ezall,2)./nanmean(xx2.eps(:,idtr),2),xx2.z,'.-','linewidth',2)
    hold on
end

xlim([0 8])
ylim([xx2.z(1) xx2.z(end)])
axis ij
grid on
%legend('0.15','0.25','0.4','0.5','0.75','1','location','best')
ylabel('Depth (m) ','fontsize',16)
xlabel('Bias Factor','fontsize',16)

%%
figdir='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/'
%fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth' ])
%save2pdfAP(fname)

%
%
%
%% Plot standard deviation of mean depth-profiles of epsilon from 1 resampling set (1 speed)

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3;
minOT=50

figure(1);clf
agutwocolumn(0.6)
wysiwyg

for testnum=[3 1  2 4 ]
    
    ezall=[];
    
    for whcase=1:100
        clear x_resamp fname
        fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
        load(fullfile('Data/',fname))
        hold on
        ezall=[ezall nanmean(x_resamp.eps,2)];
    end
    
    %~~  Load T-chain data to resample
    %     load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))
    %     % find data within time of resampling period
    %     idtr=isin(xx2.yday,[x_resamp.time(1) x_resamp.time(end)]);
    
    %    plot(nanmean(ezall,2)./nanmean(xx2.eps(:,idtr),2),xx2.z,'.-','linewidth',2)
        plot(nanstd(ezall,0,2),x_resamp.z,'.-','linewidth',2)
    
    %
    %plot(nanstd(ezall,0,2)./nanmean(ezall,2),x_resamp.z,'.-','linewidth',2)
    hold on
end

%xlim([0 8])
ylim([x_resamp.z(1) x_resamp.z(end)])
axis ij
grid on
%legend('0.15','0.25','0.4','0.5','0.75','1','location','best')
ylabel('Depth (m) ','fontsize',16)
xlabel('std','fontsize',16)

%%
figdir='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/'
%fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth' ])
%save2pdfAP(fname)
%%