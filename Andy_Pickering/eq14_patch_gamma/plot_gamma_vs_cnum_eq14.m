%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_cnum_eq14.m
%
% Plot gamma estimated for TIWE patches vs cast #
%
%
% * do for 60-200m depth also?
%
%-------------
% 2/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.4 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

eq14_patches_paths
datdir = fullfile( analysis_dir, project_long, 'data')

% load my patches
load(fullfile(datdir,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

figure(1);clf
agutwocolumn(0.75)
wysiwyg

plot(log10(patches.gam_bulk),patches.cnum,'.','color',0.85*[1 1 1])
freqline(log10(0.2))
grid on
xlabel('log_{10}\gamma_{\chi\epsilon}','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-4 1])
title(['EQ14 patches - minOT=' num2str(100*patch_size_min) 'cm' ])


% compute median gamma for each cast#
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_line(id));
end
hold on
plot(log10(gam_md),cnums,'k.')

%%

eq14_patches_paths
fname=['eq14_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gam_vs_cnum']
print(fullfile(figdir,fname),'-dpng')


%%