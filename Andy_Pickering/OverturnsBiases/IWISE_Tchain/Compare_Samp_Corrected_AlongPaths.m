%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_Samp_Corrected_AlongPaths.m
%
% Plot resampled epsilon, w/ and w/o w_t correction, along sampling paths
%
% Modifed from part of Plot_Eps_AlongSamplePaths.m
%
%---------------
% 11/23/15 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all ; clc

plotiso=1
whmoor=3 ; yl=[800 2100]
minOT=50

saveplot=1

% choose time period to plot
% tchain3
xl=[169.35 169.55];str='NearDay169' ; cl=[-8 -4]; xlprof=[1e-8 1e-5]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

diso=15

% choose which test (sampling speed)
%testnum=3 % w_samp=0.15 m/s
%testnum=2 % w_samp= 0.5 m/s
testnum=1 % w_samp=0.25 m/s

thresh=0.5 % ratio of w_t above which we call data bad

%~~~ load resampled dataset
clear fname  M N dt tc ei idsamp ep REsamp

load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

figure(1);clf
agutwocolumn(1)
wysiwyg
%ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.01, 1,3);
ax = MySubplot4(0.1, 0.03, 0.04, 0.06, 0.1, 0.03, 2,5)

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

whvar='eps' ; cl=[-8 -4]

% plot mean depth profile of true epsilon
axes(ax(1))
semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
set(gca,'Xtick',[1e-8 1e-7 1e-6 1e-5])
xlim(xlprof)
ylim(yl)
axis ij
grid on
ylabel('Depth','fontsize',16)
%SubplotLetterMW('a)')

% contour true T-chain data
axes(ax(2))
ezpc(xx2.yday(id),xx2.z,real(log10(xx2.(whvar)(:,id))))
caxis(cl)
hold on
% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
end
% plot example resampling path
plot(REsamp.tsamp(:,:,1),REsamp.zsamp(:,:,1),'w')
% make colormap for epsilon
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
xlim(xl);ylim(yl)
ytloff
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ' , w=' num2str(REsamp.w_samp) 'm/s' ],'interpreter','none')
SubplotLetterMW('true')

% color vector for plotting the resampled epsilon along paths
cvec=linspace(cl(1),cl(2),length(cmap));
Nshift=size(REsamp.eps,3);
eps_prof_all_up=nan*ones(length(REsamp.z),Nshift);
eps_prof_all_dn=nan*ones(length(REsamp.z),Nshift);


%~~ plot DOWNward sampling paths
axes(ax(4))

whdir='down'
eps_prof_all_dn=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec,whvar);
xlim(xl);ylim(yl)
ytloff
caxis([cvec(1) cvec(end)])
SubplotLetterMW('samp')

%~~ plot mean profile for downward sampling paths
axes(ax(3))
semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_dn,2),REsamp.z,'linewidth',2)
set(gca,'Xtick',[1e-8 1e-7 1e-6 1e-5])
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','sampled','location','best')
SubplotLetterMW('b)')

%~~ plot DOWNward sampling paths for CORRECTED
axes(ax(6))


whdir='down'
eps_prof_all_dn_corr=PlotEpsAlongSamplePaths_corr(REsamp,xl,whdir,cmap,cvec,thresh);

xlim(xl);ylim(yl)
ytloff
caxis([cvec(1) cvec(end)])
SubplotLetterMW('samp corr')

%~~ plot mean profile for downward sampling paths
axes(ax(5))
semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_dn,2),REsamp.z,'linewidth',2)
semilogx(nanmean(eps_prof_all_dn_corr,2),REsamp.z,'linewidth',2)
set(gca,'Xtick',[1e-8 1e-7 1e-6 1e-5])
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','samp','corr','location','best')


%~~ plot UPward sampling paths
axes(ax(8))
whdir='up'
eps_prof_all_up=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec,whvar);
xlim(xl);ylim(yl)
ytloff
caxis([cvec(1) cvec(end)])


%~~ plot mean profile for downward sampling paths
axes(ax(7))
semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_up,2),REsamp.z,'linewidth',2)
set(gca,'Xtick',[1e-8 1e-7 1e-6 1e-5])
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','sampled','location','best')


%~~ plot UPward sampling paths for CORRECTED
axes(ax(10))
whdir='up'
eps_prof_all_up_corr=PlotEpsAlongSamplePaths_corr(REsamp,xl,whdir,cmap,cvec,thresh);
xlim(xl);ylim(yl)
ytloff
caxis([cvec(1) cvec(end)])


%~~ plot mean profiles
axes(ax(9))
semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_up,2),REsamp.z,'linewidth',2)
semilogx(nanmean(eps_prof_all_up_corr,2),REsamp.z,'linewidth',2)
set(gca,'Xtick',[1e-8 1e-7 1e-6 1e-5])
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','samp','corr','location','best')
%%