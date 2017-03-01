%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_binned.m
%
%
% 3/1/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load the 1m binned Cham data processed by Sally
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam = ComputeGamma(cham.N2,cham.DTDZ_RHOORDER,cham.CHI,cham.EPSILON);

%%
figure(1);clf
histogram(log10(gam(:)),'edgecolor','none','normalization','pdf')
xlim([-4 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)
title('EQ14 1m avg \gamma_{\chi\epsilon}' )
grid on

%%
eq14_patches_paths
figdir = fullfile( analysis_dir,project,'figures')
ChkMkDir(figdir)
print(fullfile(figdir,'eq14_binned_gamma'),'-dpng')

%%