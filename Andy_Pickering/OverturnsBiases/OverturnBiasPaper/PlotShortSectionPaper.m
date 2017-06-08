%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotShortSectionPaper.m
%
% Plot a short section of Tchain data to illustrate data and some
% resampling paths etc. for the paper.
%
%
% 2 Feb 2015 - A. Pickering
% Updated June 17 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

plotiso=1
whmoor=3 % mooring #
minOT=50
saveplot=1

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps']) )


xl=[168 170]
xl=[168.3 168.6]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 183]
%xl=[180 185]
xl=[183.2 183.6]
%xl=[210 212]


figure(1);clf
set(gcf,'defaultaxesfontsize',15)
agutwocolumn(0.5)
wysiwyg
%ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 2,3);
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

ezcf(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)),cl,15)

if plotiso==1
    hold on
    %contour(xx2.yday(id),xx2.z,xx2.T(:,id),[1:7],'k')
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:12:end),'k')
end

%
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
whc=1
hold on
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w','linewidth',2)
% vline(REsamp.timeall_true(whc,:),'k--')
  vline(REsamp.tgrid(whc,:),'k--')
%  vline(REsamp.timeall(whc,:),'r--')
% testnum=5
% load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
% whc=1
% plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w','linewidth',2)
xlim(xl)

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
cmap=flipud(hot);
colormap([0.65*[1 1 1] ; cmap])
cb=colorbar
cb.Label.String='log_{10}\epsilon (Wkg^{-1}) ';
cb.FontSize=14
%
caxis(cl)
title(['Tchain ' num2str(whmoor) ])
hold off

if saveplot==1
%     figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
%     fname=fullfile(figdir,['Tchain' num2str(whmoor) '_ExampleSample' ])
%     addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
%     export_fig(fname,'-pdf')
        figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_ExampleSample' ]
    print('-dpng',fullfile(figdir,figname))

end

%%