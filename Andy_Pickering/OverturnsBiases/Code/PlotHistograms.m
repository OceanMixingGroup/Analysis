%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotHistograms.m
%
% AP June 17 2015
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot histogram of true and resampled epsilon

clear ; close all

saveplot=0
minOT=50

%~ choose mooring to plot
whmoor=3 ; xl=[1e-9 10^(-5.5)]
%whmoor=4 ; xl=[1e-9 10^(-5.9)]
%whmoor=1 ; xl=[1e-9 10^(-5.9)]
%whmoor=2 ; xl=[1e-9 10^(-5.9)]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%~ load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])

Nm='probability'
%Nm='count'
%Nm='pdf'

Nbins=40

figure(1);clf
agutwocolumn(0.5)
wysiwyg

% test # (sampling speed)
clear REsamp testnum
testnum=3
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

h1=histogram(log10(REsamp.eps_each(:)),Nbins,'Facecolor','r','Normalization',Nm)
hold on
plot(log10(nanmean((REsamp.eps_each(:)))),0,'rd','linewidth',2,'Markersize',8)
vline(log10(nanmean((REsamp.eps_each(:)))),'r--')

clear REsamp testnum
testnum=4
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

h2=histogram(log10(REsamp.eps_each(:)),h1.BinEdges,'Facecolor','b','Normalization',Nm)
plot(log10(nanmean((REsamp.eps_each(:)))),0,'bd','linewidth',2,'Markersize',8)
vline(log10(nanmean((REsamp.eps_each(:)))),'b--')

% plot `true' data
hold on
h3=histogram(log10(xx2.eps_each(:)),h1.BinEdges,'Facecolor','g','Normalization',Nm)
legend([h1 h2 h3],'0.15 m/s','0.75 m/s','true')
plot(log10(nanmean((xx2.eps_each(:)))),0,'gd','linewidth',2,'Markersize',8)
vline(log10(nanmean((xx2.eps_each(:)))),'g--')


xlim([-11 -3])
ylim([0 0.015])
ylim([0 0.1])
ylabel('probability','fontsize',16)
xlabel(['log_{10}\epsilon [Wkg^{-1}]'],'fontsize',16)
grid on
%

if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_eps_Histograms' ]
    print('-dpng',fullfile(figdir,figname))
end

%% make a 2X2 figure with histograms of Lot,Lt,N2,and eps for true and sampled quanitities

clear ; close all

saveplot=1
minOT=50
whmoor=3 ; xl=[1e-9 10^(-5.5)]
%whmoor=4 ; xl=[1e-9 10^(-5.9)]
%whmoor=1 ; xl=[1e-9 10^(-5.9)]
%whmoor=2 ; xl=[1e-9 10^(-5.9)]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
%
testnum=3
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%Nm='count'
Nm='probability'
Ds='stairs'
Ds='bar'

whvar='Lt_each'
%whvar='Lot_each'
%whvar='Otnsq_each'
%whvar='eps_each'

Nbins=20

figure(2);clf
agutwocolumn(0.8)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.05, 0.06, 0.1, 0.075, 2,2);

colsamp='m'
coltrue='g'

axes(ax(1))
whvar='Lot_each'

hs=histogram(REsamp.(whvar)(:),Nbins,'Normalization',Nm,'Facecolor',colsamp)
hold on
ht=histogram(xx2.(whvar)(:),hs.BinEdges,'Normalization',Nm,'Facecolor',coltrue)

h=vline(nanmean(REsamp.(whvar)(:)),'--')
set(h,'color',colsamp)
h=vline(nanmean(xx2.(whvar)(:)),'--')
set(h,'color',coltrue)

legend([hs ht],'samp','true','location','best')
xlabel('L')
text(600,0.3,['max=' num2str(nanmax(xx2.(whvar)(:)))],'color',coltrue,'fontweight','bold')
text(600,0.25,['max=' num2str(nanmax(REsamp.(whvar)(:)))],'color',colsamp,'fontweight','bold')

axes(ax(2))
whvar='Lt_each'

hs=histogram(REsamp.(whvar)(:),Nbins,'Normalization',Nm,'Facecolor',colsamp)
hold on
ht=histogram(xx2.(whvar)(:),hs.BinEdges,'Normalization',Nm,'Facecolor',coltrue)

h=vline(nanmean(REsamp.(whvar)(:)),'--')
set(h,'color',colsamp)
h=vline(nanmean(xx2.(whvar)(:)),'--')
set(h,'color',coltrue)

xlabel('L_T')
text(300,0.3,['max=' num2str(round(nanmax(xx2.(whvar)(:))))],'color',coltrue,'fontweight','bold')
text(300,0.25,['max=' num2str(round(nanmax(REsamp.(whvar)(:))))],'color',colsamp,'fontweight','bold')

axes(ax(3))
whvar='Otnsq_each'

hs=histogram(log10(REsamp.(whvar)(:)),Nbins,'Normalization',Nm,'Facecolor',colsamp)
hold on
histogram(log10(xx2.(whvar)(:)),hs.BinEdges,'Normalization',Nm,'Facecolor',coltrue)

h=vline(log10(nanmean(REsamp.(whvar)(:))),'--')
set(h,'color',colsamp)
h=vline(log10(nanmean(xx2.(whvar)(:))),'--')
set(h,'color',coltrue)
xlabel('log_{10}N^2')

axes(ax(4))
whvar='eps_each'

hs=histogram(log10(REsamp.(whvar)(:)),Nbins,'Normalization',Nm,'Facecolor',colsamp)
hold on
histogram(log10(xx2.(whvar)(:)),hs.BinEdges,'Normalization',Nm,'Facecolor',coltrue)
h=vline(log10(nanmean(REsamp.(whvar)(:))),'--')
set(h,'color',colsamp)
h=vline(log10(nanmean(xx2.(whvar)(:))),'--')
set(h,'color',coltrue)
xlabel('log_{10}\epsilon [Wkg^{-1}]]')

if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_2X2_Histograms' ]
    print('-dpng',fullfile(figdir,figname))
end


%%
