
%%

clear ; close all

plotiso=1
whmoor=3
minOT=50
testnum=3

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

REsamp.t=REsamp.data_resamp;
%%

em2=nanmean(REsamp.eps,3);

clear em1
for whc=1:REsamp.Nshift
    em1(:,whc)=nanmean(REsamp.eps(:,:,whc),2);
end


clear emup
for whc=1:REsamp.Nshift
    emup(:,whc)=nanmean(REsamp.eps(:,1:2:end,whc),2);
end


clear emdn
for whc=1:REsamp.Nshift
    emdn(:,whc)=nanmean(REsamp.eps(:,2:2:end,whc),2);
end

figure(1);clf
semilogx(em1,REsamp.z,'.','color',0.75*[1 1 1])
hold on
semilogx(nanmean(em1,2),REsamp.z,'k')
semilogx(nanmean(emup,2),REsamp.z)
semilogx(nanmean(emdn,2),REsamp.z)
%semilogx(nanmean(em2,2),REsamp.z,'r--')
grid on
axis ij
xlim([1e-9 1e-5])
%legend('all','up','down')

%%

figure(1);clf

for whc=1:REsamp.Nshift
    figure(1);clf
    semilogx(nanmean(REsamp.eps(:,:,whc),2),REsamp.z,'k')
    hold on
    semilogx(nanmean(REsamp.eps(:,1:2:end,whc),2),REsamp.z)
    semilogx(nanmean(REsamp.eps(:,2:2:end,whc),2),REsamp.z)
    grid on
    axis ij
    xlim([1e-9 1e-5])
    legend('all','up','down')
    pause(0.2)
end
%%