%~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotMLEcomparison.m
%
% Formerly was part of Do_MLE_Fit_EQ14cham.
%
%------------------
% 06/24/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~
%% 

clear ; close all

saveplot=0

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/

% Load data from MLE fits (see Do_MLE_Fit_EQ14cham.m)
load('/Users/Andy/Cruises_Research/ChiPod/EQ14/mfiles/MLEfits.mat')

% Load chi-pod method data
SetParamsPaperPlots

%~ load chi-pod processed data
load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));

% Bin-average Chameleon data
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)

%%

cl=[0 80]
xbins=-12:0.05:-4;
ybins=xbins;

% chipod versus Chameleon
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_bin.chi(:)),0,log10(C.chi(:)),0,0);

Nm='pdf'
%Nm='Probability'

figure(1);clf
agutwocolumn(0.8)
wysiwyg

subplot(221)
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
grid on
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
caxis(cl)
hold  on
loglog(xbins,xbins,'k--')
xlabel('\chi_{\epsilon}^{cham} [K^2s^{-1}]','fontsize',16)
ylabel('\chi_{\chi}^{cham} [K^2s^{-1}]','fontsize',16)
axis square

% Find indices for different ranges of chameleon epsilon
clear id1 id2
id1=find(log10(cham_bin.chi)<-8);
id2=find(log10(cham_bin.chi)>-8);

rat1=C.chi(id1)./cham_bin.chi(id1);
rat2=C.chi(id2)./cham_bin.chi(id2);

subplot(222)
h1=histogram(log10(rat1),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(log10(rat2),'edgecolor','none','Normalization',Nm);
freqline(nanmean(log10(rat1)),'b--')
freqline(nanmean(log10(rat2)),'r--')
xlim([-4 2])
grid on
xlabel('log_{10}[\chi_{\chi}/\chi_{\epsilon}^{cham}]','fontsize',16)
%title(['\mu_1-\mu_2=' num2str(nanmean(log10(rat1))-nanmean(log10(rat2)))])
title(['\mu_1= ' num2str(roundx(nanmean(log10(rat1)),2)) ' , \mu_2= ' num2str(roundx(nanmean(log10(rat2)),2))])
leg=legend([h1 h2],'\epsilon_1<-10e-8','\epsilon_2>10e-8','location','northwest');
leg.FontSize=14;
%
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_bin.chi(:)),0,log10(mle_bin.chi(:)),0,0);

subplot(223)
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
grid on
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
caxis(cl)
hold  on
loglog(xbins,xbins,'k--')
xlabel('\chi_{\epsilon}^{cham} [K^2s^{-1}]','fontsize',16)
ylabel('MLE method \chi [K^2s^{-1}]','fontsize',16)
axis square

clear rat1 rat2
rat1=mle_bin.chi(id1) ./ cham_bin.chi(id1);
rat2=mle_bin.chi(id2) ./ cham_bin.chi(id2);

subplot(224)
h1=histogram(log10(rat1),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(log10(rat2),'edgecolor','none','Normalization',Nm);
freqline(nanmean(log10(rat1)),'b--')
freqline(nanmean(log10(rat2)),'r--')
%leg=legend([h1 h2],'\epsilon_1<-10e-9','\epsilon_2>10e-7','location','northwest');
%leg.FontSize=14;
xlim([-4 2])
xlabel('log_{10}[\chi_{MLE}/\chi_{\epsilon}^{cham}]','fontsize',16)
grid on
title(['\mu_1= ' num2str(roundx(nanmean(log10(rat1)),2)) ' , \mu_2= ' num2str(roundx(nanmean(log10(rat2)),2))])

if saveplot==1
    figname='MLEfitsEQ14_Whist'
    SetPaperFigPath
    print( fullfile(figdir ,figname ) , '-dpng' )
end

%%
