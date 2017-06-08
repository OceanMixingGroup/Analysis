function patches=load_patches_comb(project_name,patch_size_min, usetemp, merge_patches, min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% load_patches_comb
%
%--------------
% 3/22/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
addpath(fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering',[project_name '_patch_gamma'],'code'))
eval([project_name '_patches_paths'])
%datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data';
datdir = fullfile(analysis_dir,project_long,'data');
clear patches
if merge_patches==1
    load(fullfile(datdir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) );
else
    load(fullfile(datdir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) );
end

%%