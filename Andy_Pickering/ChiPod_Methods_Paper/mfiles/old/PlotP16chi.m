%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotP16chhi.m
%
%
%
%--------
% 05/10/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/ChiPod/P16N/mfiles/

BaseDir='/Users/Andy/Cruises_Research/ChiPod/P16N/';

cruise='leg1'
load(fullfile(BaseDir,'data',['P16N_' cruise '_XC.mat']))
X1=XC.comb;clear XC

clear cruise XC P16
cruise='leg2'
load(fullfile(BaseDir,'data',['P16N_' cruise '_XC.mat']))
X2=XC.comb;clear XC

load('/Users/Andy/Cruises_Research/ChiPod/P16N/data/P16N_Locs_Depths.mat')

xl=[nanmin(P16N.leg1.lat) nanmax(P16N.leg2.lat)]
xl=[nanmin(P16N.leg2.lat) nanmax(P16N.leg2.lat)]
%xl=5*[-1 1]
%yl=[3000 5200]
%yl=[0 1000]
yl=[0 5900]
%
figure(1);clf
agutwocolumn(1)
%orient landscape
wysiwyg
ax = MySubplot(0.15, 0.075, 0.02, 0.06, 0.1, 0.02, 1,2);

set(gcf,'defaultaxesfontsize',14)

%
axes(ax(1))
%ezpc(X1.lat,X1.P,log10(X1.dTdz))
ezpc(X2.lat,X2.P,log10(X2.dTdz))
cb=colorbar
cb.Label.String='log_{10}dT/dz';
cb.Label.FontSize=16;
hold on

PlotP16BathyCombo(P16N)
xlim(xl)
ylim(yl)
xtloff
%h=SubplotLetterMW('dT/dz',1,10);;set(h,'fontsize',18)
caxis([-6 0])
plot(X2.lat,0,'k.')
title('P16N \chi pod Data')
ylabel('pressure [db]','fontsize',16)

% dtdzm=nanmean(log10(X2.dTdz),2);
% dtdzm=dtdzm(~isnan(dtdzm));
% hold on
% ig=find(~isnan(X2.lat));
% ig=ig(1:70);
% contour(X2.lat(ig),X2.P,log10(X2.dTdz(:,ig)),dtdzm(1:150:end),'k')

axes(ax(2))
%ezpc(X1.lat,X1.P,log10(X1.chi))
ezpc(X2.lat,X2.P,log10(X2.chi))
cb=colorbar;
cb.Label.String='log_{10}\chi';
cb.Label.FontSize=16;
hold on
PlotP16BathyCombo(P16N)
plot(X2.lat,0,'k.')
xlim(xl)
ylim(yl)
%xtloff
%h=SubplotLetterMW('\chi');set(h,'fontsize',18)
caxis([-11 -7])
%title('P16N Leg1 & Leg2 ')
xlabel('Latitude','fontsize',16)
ylabel('pressure [db]','fontsize',16)

linkaxes(ax)

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper'
print(fullfile(figdir,'P16N_chi'),'-dpng')

%%

figure(2);clf
ax = MySubplot(0.15, 0.075, 0.02, 0.06, 0.1, 0.02, 2,1);

axes(ax(1))
semilogx(nanmean(X2.dTdz,2),X2.P)
xlim([1e-6 1e-1])
axis ij
grid on
ylabel('Pressure [db]','fontsize',16)
xlabel('<dT/dz>','fontsize',16)

axes(ax(2))
semilogx(nanmean(smooth(X2.chi,4),2),X2.P,'.')
xlim([1e-11 1e-6])
axis ij
ytloff
grid on
xlabel('<\chi>','fontsize',16)

linkaxes(ax,'y')
%%