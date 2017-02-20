%
% Misc_Feb15.m
%
%
%
%%

clear ; close all

patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
datdir='/Users/Andy/Cruises_Research/ChiPod/TIWE/data'

load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )


%%

figure(1);clf
h1=histogram(log10(patches.dtdz_bulk(:)))
hold on
h2=histogram(log10(patches.dtdz_line(:)),h1.BinEdges)
xlabel('log_{10}[dT/dz]')

%%

figure(1);clf
h1=histogram(log10(patches.n2_bulk(:)),'Normalization','pdf');
hold on
h2=histogram(log10(patches.n2_line(:)),h1.BinEdges,'Normalization','pdf');
%h2=histogram(log10(patches.n2_ot(:)),h1.BinEdges,'Normalization','pdf');
%h3=histogram(log10(patches.n4(:)),h1.BinEdges,'Normalization','pdf');
%h2=histogram(log10(patches.n2_range(:)),h1.BinEdges,'Normalization','pdf');
xlabel('log_{10}[N^2]')
legend([h1 h2],'bulk','line')

%%

gam2=ComputeGamma(patches.n2_line,patches.dtdz_line,patches.chi_bin,patches.eps_bin);
gam3=ComputeGamma(patches.n2_bulk,patches.dtdz_bulk,patches.chi_bin,patches.eps_bin);
figure(1);clf
histogram(real(log10(patches.gam_bin)))
hold on 
histogram(log10(gam2))
histogram(log10(gam3))
freqline(log10(0.2))