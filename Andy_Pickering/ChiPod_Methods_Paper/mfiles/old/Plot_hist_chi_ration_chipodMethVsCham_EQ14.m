%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_hist_chi_ration_chipodMethVsCham_EQ14.m
%
% * Makes plot for chipod methods paper *
%
% Was part of Plot_chiProc_vs_actual_EQ14.m before
%
%---------------
% 04/06/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

saveplot=0

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures'
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/
Params.fmax=7;
Params.z_smooth=10;
Params.gamma=0.2;
Params.resp_corr=0;
Params.fc=99

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));

% load processed chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum.mat')

% average chameleon data in same depth bins as chipod method
% so we can make scatter plots etc further on
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)


%% Plot histogram of chipod chi / cham chj

C.chi(find(log10(C.chi)>-2))=nan;
C.chi(find(log10(C.chi)<-13))=nan;

cham_bin.chi(find(log10(cham_bin.chi)>-2))=nan;
cham_bin.chi(find(log10(cham_bin.chi)<-13))=nan;

figure(1);clf
chirat=C.chi(:)./cham_bin.chi(:);
histogram(log10(chirat))
xlim(2*[-1 1])
xlabel('log_{10} [\chi_{\chi}/\chi_{\epsilon}]','fontsize',16)
grid on
hold on

line(nanmean(log10(chirat))*[1 1],[0 5000],'color','k')
plot(nanmean(log10(chirat)),0,'kd','linewidth',2)

line(nanmedian(log10(chirat))*[1 1],[0 5000],'color','r')
plot(nanmedian(log10(chirat)),0,'rd','linewidth',2)

ylabel('count','fontsize',16)

if saveplot==1
    figname=['EQ14_chi_ratio_hist_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
    print( fullfile('/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper' ,figname ) , '-dpng' )
end

%%