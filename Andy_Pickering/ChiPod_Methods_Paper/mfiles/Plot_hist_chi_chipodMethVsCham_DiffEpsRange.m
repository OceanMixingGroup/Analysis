%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_hist_chi_chipodMethVsCham_DiffEpsRange.m
%
% Look at the ratio of chi_chipod / chi:cham for different ranges of
% epsilon (as measured by chameleon).
%
% see also Plot_2Dhist_chi_chipodMethVsCham.m
%
%---------------
% 08/03/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

saveplot=0

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

% Set chipod params and load chameleon data
SetParamsPaperPlots
Params

%~ load chi-pod processed data
load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));
%
% average chameleon data in same depth bins as chipod method
% so we can make scatter plots etc further on
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)

% Same as above, with histogram added

saveplot=0

clear xbins ybins
%xbins=-12:0.1:-4;
xbins=-12:0.05:-4;
ybins=xbins;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_bin.chi(:)),0,log10(C.chi(:)),0,0);
%
figure(9);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

%subplot(1,3,[1 2])
subplot(211)
h=pcolor(xbins,ybins,hist);
set(h,'edgecolor','none')
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
%axis square
grid on
cmap=flipud(hot);
colormap([0.85*[1 1 1] ; cmap])
caxis([0 80])
ylabel('log_{10}\chi_{\chi}^{cham} [K^2s^{-1}]','fontsize',16)
xlabel('log_{10}\chi_{\epsilon}^{cham} [K^2s^{-1}]','fontsize',16)
% ylabel('\chi pod method : log_{10}\chi_{\chi}^{cham} [K^2s^{-1}]','fontsize',16)
% xlabel('Chameleon : log_{10}\chi_{\epsilon}^{cham} [K^2s^{-1}]','fontsize',16)
title('2D hist of \chi from both methods','fontsize',16)
set(gca,'Xtick',[-12:1:-4])

%subplot(1,3,3)
subplot(212)

chirat=C.chi(:) ./  cham_bin.chi(:);
histogram(log10( chirat ) ,'EdgeColor','none','normalization','pdf','facecolor',[254 69 0]./255)
xlim(2.2*[-1 1])
grid on
vline(log10(nanmedian(chirat)),'k--')
xlabel('log_{10}[\chi_{\chi}^{cham} / \chi_{\epsilon}^{cham}]','fontsize',16)
title(['\mu=' num2str(roundx(nanmean(log10(chirat)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat)),2))])
set(gca,'Xtick',[-2:0.5:2])
ylabel('counts','fontsize',16)

% if saveplot==1
%     figname=['EQ14_2dhist_chi_WITHhist_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
%     SetPaperFigPath
%     print( fullfile(figdir ,figname ) , '-dpng' )
% end

%% make histograms for different ranges of epsilon

clear id id1 id2 

id=find( log10(cham_bin.eps)<-10);
chirat=C.chi(id) ./  cham_bin.chi(id);

id1=find( log10(cham_bin.eps)>-10 & log10(cham_bin.eps)<-8);
chirat1=C.chi(id1) ./  cham_bin.chi(id1);

id2=find( log10(cham_bin.eps)>-8 & log10(cham_bin.eps)<-6);
chirat2=C.chi(id2) ./  cham_bin.chi(id2);

id3=find( log10(cham_bin.eps)>-6);
chirat3=C.chi(id3) ./  cham_bin.chi(id3);

rr=3
cc=1

figure(2);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(rr,cc,1)
h1=histogram(log10( chirat1 ),'BinEdges',[-2:0.05:2] ,'EdgeColor','none','normalization','pdf','facecolor',[254 69 0]./255)
hold on
vline(log10(nanmedian(chirat1)),'k--')
xlim(2.2*[-1 1])
grid on
title(['-10 < log_{10}\epsilon < -8, \mu=' num2str(roundx(nanmean(log10(chirat1)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat1)),2))])

ax2=subplot(rr,cc,2)
histogram(log10( chirat2 ),h1.BinEdges ,'EdgeColor','none','normalization','pdf')
hold on
vline(log10(nanmedian(chirat2)),'k--')
xlim(2.2*[-1 1])
grid on
title(['-8 < log_{10}\epsilon < -6, \mu=' num2str(roundx(nanmean(log10(chirat2)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat2)),2))])

ax3=subplot(rr,cc,3)
histogram(log10( chirat3 ),h1.BinEdges ,'EdgeColor','none','normalization','pdf')
hold on
vline(log10(nanmedian(chirat3)),'k--')
xlim(2.2*[-1 1])
grid on
title(['-6 < log_{10}\epsilon , \mu=' num2str(roundx(nanmean(log10(chirat3)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat3)),2))])

% ax4=subplot(rr,cc,4)
% histogram(log10( chirat3 ),h1.BinEdges ,'EdgeColor','none','normalization','pdf')
% hold on
% vline(log10(nanmedian(chirat3)),'k--')
% xlim(2.2*[-1 1])
% grid on

xlabel('log_{10}[\chi_{\chi}^{cham} / \chi_{\epsilon}^{cham}]','fontsize',16)


linkaxes([ax1 ax2 ax3 ])
%%

vline(log10(nanmedian(chirat)),'k--')
xlabel('log_{10}[\chi_{\chi}^{cham} / \chi_{\epsilon}^{cham}]','fontsize',16)
title(['\mu=' num2str(roundx(nanmean(log10(chirat)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat)),2))])
set(gca,'Xtick',[-2:0.5:2])
ylabel('counts','fontsize',16)



%%
% %% 2D hist of chi from both methods
% 
% clear xbins ybins
% xbins=-12:0.05:-4;
% ybins=xbins;
% 
% addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
% [hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_bin.chi(:)),0,log10(C.chi(:)),0,0);
% %
% figure(9);clf
% agutwocolumn(0.5)
% wysiwyg
% h=pcolor(xbins,ybins,hist);
% set(h,'edgecolor','none')
% hold on
% plot(xbins,xbins,'k--')
% plot(xbins,xbins-1,'r--')
% plot(xbins,xbins+1,'r--')
% colorbar
% grid on
% cmap=flipud(hot);
% colormap([0.85*[1 1 1] ; cmap])
% caxis([0 80])
% ylabel('\chi pod method : log_{10}\chi_{\chi} ','fontsize',16)
% xlabel('chameleon : log_{10}\chi_{\epsilon} ','fontsize',16)
% title('2D hist of \chi from both methods','fontsize',16)
% 
% if saveplot==1
%     figname=['EQ14_2dhist_chi_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
%     SetPaperFigPath
%     print( fullfile(figdir ,figname ) , '-dpng' )
% end

%%
