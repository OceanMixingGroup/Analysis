%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_depth_eq14.m
%
%
%-----------
% 3/1/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam = ComputeGamma(cham.N2,cham.DTDZ_RHOORDER,cham.CHI,cham.EPSILON);

%%

figure(1);clf
agutwocolumn(0.6)
wysiwyg
%h1=plot(log10(gam(:)),cham.P(:),'.')
histogram2(log10(gam(:)),cham.P(:),'displaystyle','tile')
axis ij
ylim([0 250])
grid on
xlim([-5 2])
ylabel('P')
xlabel('log_{10}[\gamma]')
freqline(log10(0.2))
cmap=colormap;
colormap([0.75*[1 1 1];cmap])
title('EQ14 binned \gamma vs depth')

%% Do same for patch gamma

clear ; close all

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')

figure(1);clf
%plot(log10(patches.gam_line(:)),patches.p1(:),'.')
histogram2(log10(patches.gam_line(:)),patches.p1(:),'displaystyle','tile')
axis ij
ylim([0 250])
grid on
xlim([-5 2])
ylabel('P')
xlabel('log_{10}[\gamma]')
freqline(log10(0.2))
cmap=colormap;
colormap([0.75*[1 1 1];cmap])
title('EQ14 patch \gamma vs depth')
caxis([0 1000])


%% Plot depth-distribution of patches

figure(2);clf
histogram(nanmean([patches.p1(:) patches.p2(:)],2))
xlabel('patch depth (mean)')
ylabel('count')
title('EQ14 patch depth distribution')

%%