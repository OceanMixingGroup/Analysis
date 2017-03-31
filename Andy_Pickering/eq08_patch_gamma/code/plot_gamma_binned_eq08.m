%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_binned_eq08.m
%
% Plot histogram of the mixing coefficient estimated from the 1m avg
% Chameleon data for EQ08
%
%
%----------------
% 3/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load the 1m binned Cham data
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
cham1 = cham; clear cham

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam1 = ComputeGamma(cham1.N2,cham1.DTDZ,cham1.CHI,cham1.EPSILON);
%%
%z_range=[0 200]
z_range = [60 200]
id=find(cham1.P>z_range(1) & cham1.P<z_range(2));
%%

figure(1);clf
histogram(log10(gam1(id)),'edgecolor','none','normalization','pdf')
hold on
%histogram(log10(gam2(:)),'edgecolor','none','normalization','pdf')
xlim([-4 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)
title('EQ08 1m avg \gamma_{\chi\epsilon}' )
grid on

%%

eq08_patches_paths
print(fullfile(fig_dir,'eq08_binned_gamma'),'-dpng')

%%
