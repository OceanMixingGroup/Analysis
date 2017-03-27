%%
%
% Plot_eps_vs_chi.m
%
% for EQ14?
%
%
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

figure(1);clf
agutwocolumn(0.75)
wysiwyg

subplot(121)
histogram2(log10(cham.CHI(:)), log10(cham.EPSILON(:)),'DisplayStyle','tile')
xlim([-12 -4])
ylim([-12 -4])
xlabel('log_{10}[\chi]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)
title('all depths')

ig = find(cham.P>80);
subplot(122)
histogram2(log10(cham.CHI(ig)), log10(cham.EPSILON(ig)),'DisplayStyle','tile')
xlim([-12 -4])
ylim([-12 -4])
xlabel('log_{10}[\chi]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)

title('below 80m')

eq14_patches_paths
print(fullfile(fig_dir,'eq14_epsVschi_1mbinned'),'-dpng')

%%


%%