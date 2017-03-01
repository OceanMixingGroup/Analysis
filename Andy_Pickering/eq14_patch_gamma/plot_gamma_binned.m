%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_binned.m
%
% Plot histogram of the mixing coefficient estimated from the 1m avg
% Chameleon data for EQ14
%
% 3/1/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load the 1m binned Cham data processed by Sally
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
cham1 = cham; clear cham

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')
cham2 = cham; clear cham

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam1 = ComputeGamma(cham1.N2,cham1.DTDZ_RHOORDER,cham1.CHI,cham1.EPSILON);
%gam2 = ComputeGamma(cham2.N2,cham2.DTDZ_RHOORDER,cham2.CHI,cham2.EPSILON);
%%
figure(1);clf
histogram(log10(gam1(:)),'edgecolor','none','normalization','pdf')
hold on
%histogram(log10(gam2(:)),'edgecolor','none','normalization','pdf')
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

%% Plot estimates for fmax=7hz ?

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')
cham1=cham; clear cham

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')
cham2 = cham; clear cham

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam1 = ComputeGamma(cham1.N2,cham1.DTDZ_RHOORDER,cham1.CHI,cham1.EPSILON);
gam2 = ComputeGamma(cham2.N2,cham2.DTDZ_RHOORDER,cham2.CHI,cham2.EPSILON);

figure(1);clf
h1=histogram(log10(gam1(:)),'edgecolor','none','normalization','pdf')
hold on
h2=histogram(log10(gam2(:)),'edgecolor','none','normalization','pdf')
xlim([-4 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)
title('EQ14 1m avg \gamma_{\chi\epsilon}' )
grid on
legend([h1 h2],'7hz','32hz')

%%
eq14_patches_paths
figdir = fullfile( analysis_dir,project,'figures')
ChkMkDir(figdir)
print(fullfile(figdir,'eq14_binned_gamma_7hz'),'-dpng')

%%

figure(2);clf
histogram(log10(gam1(:)./gam2(:)),'edgecolor','none')
xlim([-1 1])

%%
