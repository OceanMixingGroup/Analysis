%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotEpsProfiles2X2_v2.m
%
%
%
%~~~~~~~~~~~~~~
% June 17, 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=0
minOT=50
whmoor=3 ; xl=[1e-9 10^(-5.5)]
%whmoor=4 ; xl=[1e-9 10^(-5.9)]
%whmoor=1 ; xl=[1e-9 10^(-5.9)]
%whmoor=2 ; xl=[1e-9 10^(-5.9)]

cd /Users/Andy/Cruises_Research/AnalysisAndy_Pickering/OverturnsBiases

% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
%

ezall1=[];
testnum=1
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall1=squeeze(nanmean(REsamp.eps,2));

%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);
zvec=REsamp.z;
clear REsamp

ezall2=[];
testnum=2
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall2=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall3=[];
testnum=3
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall3=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall4=[];
testnum=4
% load resampled data
clear REsamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall4=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall5=[];
testnum=5
% load resampled data
clear REsamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall5=squeeze(nanmean(REsamp.eps,2));
clear REsamp

%%
myc=[0.8500    0.3250    0.0980];
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 2,2);

axes(ax(1))
% plot resampled mean profiles

interval=100
nboot=100
addpath /Users/Andy/Cruises_Research/ChiPod/IWISE10_And_vmp/
[xout bb hm]=BinAndBootProfiles(ezall3,zvec,interval,nboot,1);
hold on
% plot actual mean profile
[xout bb hm]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);

% plot +/- standard devation
%semilogx(nanmean(ezall3,2)-std(ezall3,0,2),zvec,'--','color',0.3*[1 1 1])
%semilogx(nanmean(ezall3,2)+std(ezall3,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
axis ij
ylabel('Depth ','fontsize',16)
ht=SubplotLetter(['w=0.15 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on

%

axes(ax(2))

[xout bb hm]=BinAndBootProfiles(ezall1,zvec,interval,nboot,1);
hold on
% plot actual mean profile
[xout bb hm]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);


% plot resampled mean profiles
% semilogx(ezall1,zvec,'.','color',0.7*[1 1 1])
% hold on
% h1=semilogx(nanmean(ezall1,2),zvec,'linewidth',2,'color',myc)
% hold on
% % plot actual mean profile
% ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% % plot +/- standard devation
% semilogx(nanmean(ezall1,2)-std(ezall1,0,2),zvec,'--','color',0.3*[1 1 1])
% semilogx(nanmean(ezall1,2)+std(ezall1,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
ytloff
axis ij
ht=SubplotLetter(['w=0.25 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on

axes(ax(3))
[xout bb hm]=BinAndBootProfiles(ezall2,zvec,interval,nboot,1);
hold on
% plot actual mean profile
[xout bb hm]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);

% semilogx(ezall2,zvec,'.','color',0.7*[1 1 1])
% hold on
% % plot resampled mean profiles
% h2=semilogx(nanmean(ezall2,2),zvec,'linewidth',2,'color',myc)
% hold on
% % plot actual mean profile
% ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% % plot +/- standard devation
% semilogx(nanmean(ezall2,2)-std(ezall2,0,2),zvec,'--','color',0.3*[1 1 1])
% semilogx(nanmean(ezall2,2)+std(ezall2,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
axis ij
xlabel('<\epsilon> (Wkg^{-1})','fontsize',16)
ylabel('Depth ','fontsize',16)
ht=SubplotLetter(['w=0.5 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on
%
axes(ax(4))
[xout bb h1]=BinAndBootProfiles(ezall4,zvec,interval,nboot,1);
hold on
% plot actual mean profile
[xout bb h2]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);
legend([h1 h2],'samp','true')

% semilogx(ezall4,zvec,'.','color',0.7*[1 1 1])
% hold on
% % plot resampled mean profiles
% h4=semilogx(nanmean(ezall4,2),zvec,'linewidth',2,'color',myc)
% hold on
% % plot actual mean profile
% ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% % plot +/- standard devation
% semilogx(nanmean(ezall4,2)-std(ezall4,0,2),zvec,'--','color',0.3*[1 1 1])
% semilogx(nanmean(ezall4,2)+std(ezall4,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
ytloff
axis ij
xlabel('<\epsilon> (Wkg^{-1})','fontsize',16)
ht=SubplotLetter(['w=0.75 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on
shg

linkaxes(ax)

%%

figure(2);clf

interval=50
nboot=100
addpath /Users/Andy/Cruises_Research/ChiPod/IWISE10_And_vmp/

[xout bb h1]=BinAndBootProfiles(ezall3,zvec,interval,nboot,1);
hold on

[xout bb h2]=BinAndBootProfiles(ezall1,zvec,interval,nboot,1);
[xout bb h3]=BinAndBootProfiles(ezall2,zvec,interval,nboot,1);
%[xout bb h4]=BinAndBootProfiles(ezall4,zvec,interval,nboot,1);
[xout bb h4]=BinAndBootProfiles(ezall5,zvec,interval,nboot,1);
% plot actual mean profile
[xout bb h5]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);

legend([h1 h2 h3 h4 h5])

%%