%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% plot_patch_locations_tiwe.m
%
%
%-----------
% 3/2/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplots=1

% patch options
patch_size_min = 0.15 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

% option to use merged patches
merge_patches = 1 ;
min_sep = 0.15 ;


% load binned Chameleon data (to plot epsilon etc.)
load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_1mavg_combined.mat')

figure(1);clf
agutwocolumn(0.75)
wysiwyg

ezpc(cham.cnum,cham.P,log10(cham.EPSILON))
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon'
ylabel('P')
xlabel('cast #')


% load my patches
if merge_patches==1
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) )
else
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) )
end

hold on
plot(patches.cnum,patches.p1,'k.')
%plot(patches.cnum,patches.p2,'.')
xlim([2800 3700])
xlim([3000 3400])

%%

if saveplots==1
fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures';

if merge_patches==1
    fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patch_locs_merged_minsep_' num2str(min_sep*100)]
else
    fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patch_locs' ]
end

print(fullfile(fig_dir,fname),'-dpng')
end


%%
