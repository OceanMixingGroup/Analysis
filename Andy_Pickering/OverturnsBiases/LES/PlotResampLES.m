%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotResampLES.m
%
% Quick and dirty plot of mean profiles from true and resampled overturns
%
% 20 Feb 2015
% 24 Feb
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

testnum=3
minOT=50

Bdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/Data'
load(fullfile(Bdir,['Test' num2str(testnum)],['LES_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases.mat']))

figure(1);clf

Nshift=size(REsamp.eps,3)
emp=nan*ones(length(REsamp.z),Nshift);

for whc=1:Nshift
   semilogx(nanmean(REsamp.eps(:,:,whc),2),REsamp.z,'.','color',0.7*[1 1 1])
   hold on
   emp(:,whc)=nanmean(REsamp.eps(:,:,whc),2);
end

axis ij
grid on
hs=semilogx(nanmean(emp,2),REsamp.z,'k','linewidth',2)
xlim([1e-9 10^(-5.5)])
ylabel('Depth (m) ')
xlabel('\epsilon')
ylim([nanmin(REsamp.z) nanmax(REsamp.z)])
%
% plot 'true' profile
load(['N2_OT_minOT_' num2str(minOT) '.mat'])
ht=semilogx(nanmean(OT.eps,2),OT.z,'r')
shg

legend([ht hs],'true','samp','location','best')

%% Compare depth-time plots
close all

for whc=1:size(REsamp.eps,3)

cl=[-9 -4]

figure(2);clf
agutwocolumn(0.7)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 1,2);
sgm=nanmean(OT.sgth,2);
axes(ax(1))
ezpc(OT.yday,OT.z,log10(OT.eps))
hold on
contour(OT.yday,OT.z,OT.sgth,sgm(1:25:end),'k')
caxis(cl)
colorbar
colormap(jet)
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')

axes(ax(2))
ezpc(REsamp.tgrid(whc,:),REsamp.z,log10(REsamp.eps(:,:,whc)))
hold on
contour(REsamp.tgrid(whc,:),REsamp.z,REsamp.data_resamp(:,:,whc),sgm(1:25:end),'k')
caxis(cl)
colorbar

linkaxes(ax)

pause(0.2)

end

%%