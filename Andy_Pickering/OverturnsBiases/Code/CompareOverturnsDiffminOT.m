%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareOverturnsDiffminOT.m
%
% Compare overturn quanitities using different minimum overturn size
% criteria for T-chains. Data is gridded to 10 m, but actual vertical
% spacing was more like ~50m so we probably can't resolve overturns < 50m.
%
%-----------------------
% 27 Jan 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

whmoor=4

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases


load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))
X1=xx2;clear xx2

minOT=50
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]))
X2=xx2;clear xx2

minOT=100
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]))
X3=xx2;clear xx2

minOT=250
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]))
X4=xx2;clear xx2

%% Check that histograms look same for L>minOT

Nm='probability'
%Nm='count'
figure(1);clf
histogram(X1.Lot,'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(X2.Lot,'DisplayStyle','Stair','Normalization',Nm)

%%
figure(1);clf
h1=semilogx(nanmean(X1.eps,2),X1.z)
hold on
h2=semilogx(nanmean(X2.eps,2),X1.z,'--')
h3=semilogx(nanmean(X3.eps,2),X1.z,'--')
h4=semilogx(nanmean(X4.eps,2),X1.z,'--')
axis ij
xlim([1e-9 1e-6])
legend('10m',[num2str(X2.minOT) 'm'],[num2str(X3.minOT) 'm'],[num2str(X4.minOT) 'm'])

%%