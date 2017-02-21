%~~~~~~~~~~~~~~~~~~~
%
% compare_patches_tiwe_AP_Bill.m
%
% Compare my patch estimates for TIWE to the results of Bill's analysis. I
% don't have the code he used to produce it, just the resulting patches.
%
%
%------------
% 2/21/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

% load my patches
load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

% load Bills patches
load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/events_TIWE.mat')

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
histogram(log10(patches.n2_line(:)),'Normalization','pdf')
hold on
histogram(log10(A.N2(:)),'Normalization','pdf')
xlabel('log_{10}[N^2]')

subplot(222)
histogram(log10(patches.dtdz_line(:)),'Normalization','pdf')
hold on
histogram(real(log10(A.tgrad(:))),'Normalization','pdf')
xlabel('log_{10}[T_z]')

subplot(223)
histogram(log10(patches.chi(:)),'Normalization','pdf')
hold on
histogram(real(log10(A.chi(:))),'Normalization','pdf')
xlabel('log_{10}[\chi]')

subplot(224)
histogram(log10(patches.eps(:)),'Normalization','pdf')
hold on
histogram(real(log10(A.eps(:))),'Normalization','pdf')
xlabel('log_{10}[\epsilon]')

%% compute gamma from Bill's patches

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

gam_bill = ComputeGamma(A.N2,A.tgrad,A.chi,A.eps);

figure(2);clf
histogram(log10(gam_bill(:)))
hold on
histogram(log10(patches.gam_line(:)))
freqline(log10(0.2))

%%