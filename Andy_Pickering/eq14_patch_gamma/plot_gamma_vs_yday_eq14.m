%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_yday_eq14.m
%
% Plot gamma estimated for TIWE patches vs cast #
%
%-------------
% 2/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

eq14_patches_paths
datdir = fullfile( analysis_dir, project, 'data')
%datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

% load my patches
load(fullfile(datdir,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )


figure(1);clf
agutwocolumn(0.75)
wysiwyg

plot(log10(patches.gam_bulk),patches.cnum,'.')
freqline(log10(0.2))
grid on
xlabel('log_{10}\Gamma','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-4 1])
title(['EQ14 patches - minOT=' num2str(100*patch_size_min) 'cm' ])
%
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_line(id));
end
hold on
plot(log10(gam_md),cnums,'r.')

%%
fig_dir = fullfile( analysis_dir, project, 'figures')
ChkMkDir(fig_dir)
fname=['eq14_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gam_vs_cnum']
print(fullfile(fig_dir,fname),'-dpng')


%%