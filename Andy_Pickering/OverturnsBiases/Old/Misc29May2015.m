%%



%%

clear ; close all


clear ; close all

plotiso=1
whmoor=3 ; yl=[800 2100]
%whmoor=4 ; yl=[900 2200]
minOT=50

saveplot=0

% tchain3
%xl=[169.35 169.55];str='NearDay169' ; cl=[-8 -4]; xlprof=[1e-8 1e-5]
%xl=[168.25 168.7];str='NearDay168' ; cl=[-8 -4.2]; xlprof=[1e-8 1e-5]
%xl=[182.2 182.5]; str='NearDay182' ; cl=[-8 -3.6]; xlprof=[1e-8 1e-4]
%xl=[183.3 183.5] ; str='NearDay183' ; cl=[-8 -3.6]; xlprof=[1e-8 1e-4]
xl=[198.3 198.6] ; str='NearDay198' ; cl=[-8 -4.5]; xlprof=[1e-8 1e-5]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

diso=15

%for
testnum=1%[2 3]

%~~~ load resampled dataset
clear fname  M N dt tc ei idsamp ep REsamp

load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%
figure(1);clf
agutwocolumn(1)
wysiwyg
%ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.01, 1,3);
ax = MySubplot4(0.1, 0.03, 0.04, 0.06, 0.1, 0.03, 2,3)

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

% plot true T-chain data
axes(ax(2))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
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
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ' , w=' num2str(REsamp.w_samp) 'm/s' ],'interpreter','none')
SubplotLetterMW('\epsilon true')

% plot mean depth profile of true epsilon
axes(ax(1))
semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
axis ij
grid on
ylabel('Depth','fontsize',16)

% color vector for plotting the resampled epsilon along paths
cvec=linspace(cl(1),cl(2),length(cmap));
%

axes(ax(4))

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
end

whdir='down'
%whdir='up'

eps_prof_all_dn=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec);

xlim(xl);ylim(yl)
caxis([cvec(1) cvec(end)])

% plot mean profile for downward sampling paths
axes(ax(3))
semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_dn,2),REsamp.z)
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','sampled','location','best')

%

axes(ax(6))

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
end

%whdir='down'
eps_prof_all_dn=PlotEpsAlongSamplePaths_corr(REsamp,xl,whdir,cmap,cvec);

xlim(xl);ylim(yl)
caxis([cvec(1) cvec(end)])

% plot mean profile for downward sampling paths
axes(ax(5))
semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
xlim(xlprof)
ylim(yl)
hold on
semilogx(nanmean(eps_prof_all_dn,2),REsamp.z)
axis ij
grid on
ylabel('Depth','fontsize',16)
legend('true','sampled','location','best')

%%
