
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3
minOT=50
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

%%

ig=find(log10(xx2.eps)>-10);

figure(1);clf
histogram(log10(xx2.eps(ig)))
xlim([-9 -4])
%%
Ds='stair'
Nm='probability'
%Nm='count'
Nm='pdf'

Nbins=30
idt=isin(xx2.yday,[165 185]);
figure(1);clf
histogram((xx2.Lot(:,idt)),Nbins,'DisplayStyle',Ds,'Normalization',Nm)
hold on

clear REsamp
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

histogram(REsamp.Lot(:),Nbins,'DisplayStyle',Ds,'Normalization',Nm)

clear REsamp
testnum=2
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

histogram(REsamp.Lot(:),Nbins,'DisplayStyle',Ds,'Normalization',Nm)

clear REsamp
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

histogram(REsamp.Lot(:),Nbins,'DisplayStyle',Ds,'Normalization',Nm)

%xlim([-9 -4])
%%
Ds='stair'
Nm='probability'
%Nm='count'
Nm='pdf'

%Nbins=30
idt=isin(xx2.yday,[165 185]);

eps=xx2.eps(:,idt);
ig=find(log10(eps)>-10);

figure(1);clf
agutwocolumn(0.6)
wysiwyg
histogram(log10(eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','k')
hold on
%h2=vline(nanmedian(log10(eps(ig))),'--')
%h2=vline(nanmean(log10(eps(ig))),'--')
h2=vline(log10(nanmean(eps(ig))),'--')
h2.Color='k';
h2=vline(log10(nanmedian(eps(ig))),'-.')
h2.Color='k';

clear REsamp
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
clear ig
ig=find(log10(REsamp.eps)>-10);
histogram(log10(REsamp.eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','r')
%h2=vline(nanmedian(log10(REsamp.eps(ig))),'--')
%h2=vline(nanmean(log10(REsamp.eps(ig))),'--')
h2=vline(log10(nanmean(REsamp.eps(ig))),'--')
h2.Color='r'
h2=vline(log10(nanmedian(REsamp.eps(ig))),'-.')
h2.Color='r'
% 
% clear REsamp
% testnum=2
% load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
% clear ig
% ig=find(log10(REsamp.eps)>-10);
% histogram(log10(REsamp.eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','b')
% %h2=vline(nanmedian(log10(REsamp.eps(ig))),'--')
% h2=vline(nanmean(log10(REsamp.eps(ig))),'--')
% set(h2,'color','b')
% 
% clear REsamp
% testnum=4
% load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
% 
% clear ig
% ig=find(log10(REsamp.eps)>-10);
% h=histogram(log10(REsamp.eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','m')
% %h2=vline(nanmedian(log10(REsamp.eps(ig))),'--')
% h2=vline(nanmean(log10(REsamp.eps(ig))),'--')
% set(h2,'color','m')

legend('true','0.25')

grid on
%xlim([-9 -4])
xlabel('log_{10} \epsilon')
title(['Tchain ' num2str(whmoor)])
%%

Ds='stair'
Nm='probability'
%Nm='count'
Nm='pdf'

%Nbins=30
idt=isin(xx2.yday,[165 185]);

eps=xx2.eps(:,idt);
Lot=xx2.Lot(:,idt);
ig=find(log10(eps)>-10);

figure(1);clf
agutwocolumn(0.6)
wysiwyg
histogram((Lot(ig)),30,'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','k')
hold on
h2=vline(nanmean(Lot(ig)),'--')
h2.Color='k';
h2=vline(nanmedian(Lot(ig)),'-.')
h2.Color='k';
text(400,0.006,['max=' num2str(nanmax(Lot(ig))) ],'color','k')
text(400,0.005,['mean=' num2str(nanmean(Lot(ig))) ],'color','k')
text(400,0.004,['med=' num2str(nanmedian(Lot(ig))) ],'color','k')
%
clear REsamp
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
clear ig
ig=find(log10(REsamp.eps)>-10);
histogram(REsamp.Lot(ig),30,'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','r')
h2=vline(nanmean(REsamp.Lot(ig)),'--')
h2.Color='r'
h2=vline(nanmedian(REsamp.Lot(ig)),'-.')
h2.Color='r'
nanmax(REsamp.Lot(ig))
text(600,0.006,['max=' num2str(nanmax(REsamp.Lot(ig))) ],'color','r')
text(600,0.005,['mean=' num2str(nanmean(REsamp.Lot(ig))) ],'color','r')
text(600,0.004,['med=' num2str(nanmedian(REsamp.Lot(ig))) ],'color','r')
% 
% clear REsamp
% testnum=2
% load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
% clear ig
% ig=find(log10(REsamp.eps)>-10);
% histogram(log10(REsamp.eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','b')
% %h2=vline(nanmedian(log10(REsamp.eps(ig))),'--')
% h2=vline(nanmean(log10(REsamp.eps(ig))),'--')
% set(h2,'color','b')
% 
% clear REsamp
% testnum=4
% load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
% 
% clear ig
% ig=find(log10(REsamp.eps)>-10);
% h=histogram(log10(REsamp.eps(ig)),'DisplayStyle',Ds,'Normalization',Nm,'EdgeColor','m')
% %h2=vline(nanmedian(log10(REsamp.eps(ig))),'--')
% h2=vline(nanmean(log10(REsamp.eps(ig))),'--')
% set(h2,'color','m')

legend('true','0.25')

grid on
%xlim([-9 -4])
xlabel('L (m)')
title(['Tchain ' num2str(whmoor)])


%% Plot cumulative epsilon as function of increasing Lot

clear idt Lvec B I epsvec eps_sort CS Nvec 
idt=isin(xx2.yday,[180 185]);
Lvec=xx2.Lot(:,idt);Lvec=Lvec(:);
[B,I] = sort(Lvec,1,'ascend');
%[B,I] = sort(Lvec,1,'descend');
%
epsvec=xx2.eps(:,idt);epsvec=epsvec(:);
eps_sort=epsvec(I);

CS=nancumsum(epsvec(I));
Nvec=1:numel(epsvec);
truemean=nanmean(epsvec);

%
figure(1);clf
agutwocolumn(0.6)
wysiwyg
plot(B,CS./nansum(epsvec),'-')
%plot(B,CS./numel(epsvec),'.');hline(nanmean(epsvec))
%plot(B,CS./Nvec' ,'.')
xlabel('L')
ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon')
grid on
%
hold on

clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Lvec=REsamp.Lot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));
%Nvec=1:numel(epsvec);
%truemean=nanmean(epsvec);
plot(B,CS./nansum(epsvec),'-')

clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Lvec=REsamp.Lot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));
plot(B,CS./nansum(epsvec),'-')

xlabel('L')
ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon')
grid on
title(['T-chain ' num2str(whmoor)])
legend('true','0.25m/s','0.75m/s','location','best')

%%
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_EpsContribVsLot'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')

%%