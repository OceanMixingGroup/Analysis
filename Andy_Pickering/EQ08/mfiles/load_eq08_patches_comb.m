function patches=load_eq08_patches_comb(patch_size_min, usetemp, merge_patches, min_sep)

datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq08_patch_gamma/data';
eq08_patches_paths
clear patches
if merge_patches==1
    load(fullfile(datdir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) );
else
    load(fullfile(datdir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) );
end

%%