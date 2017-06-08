%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_2Dhist_chi_chipodMethVsCham.m
%
% * Makes plot for chipod methods paper *
%
% Plot a 2D histogram of \chi from chipod method versus \chi from Chameleon
% for EQ14 Chameleon profiles.
%
% Was part of Plot_chiProc_vs_actual_EQ14.m before
%
%
%---------------
% 04/06/16 - A.Pickering - apickering@coas.oregonstate.edu
% 06/03/16 - AP - Clean up and comment.
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

% average chameleon data in same depth bins as chipod method
% so we can make scatter plots etc further on
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)

% Same as above, with histogram added

saveplot=1

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
title('2D hist of \chi from both methods','fontsize',16)
set(gca,'Xtick',[-12:1:-4])

%subplot(1,3,3)
subplot(212)

chirat=C.chi(:) ./ cham_bin.chi(:);
histogram(log10( chirat ) ,'EdgeColor','none','normalization','pdf','facecolor',[254 69 0]./255)
xlim(2.2*[-1 1])
grid on
vline(log10(nanmedian(chirat)),'k--')
xlabel('log_{10}[\chi_{\chi}^{cham} / \chi_{\epsilon}^{cham}]','fontsize',16)
title(['\mu=' num2str(roundx(nanmean(log10(chirat)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(chirat)),2))])
set(gca,'Xtick',[-2:0.5:2])
ylabel('counts','fontsize',16)

if saveplot==1
    figname=['EQ14_2dhist_chi_WITHhist_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
    SetPaperFigPath
    print( fullfile(figdir ,figname ) , '-dpng' )
end

%%
