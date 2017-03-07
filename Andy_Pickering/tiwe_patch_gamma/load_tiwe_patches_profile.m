function patches=load_tiwe_patches_profile(cnum, patch_size_min,...
    usetemp, merge_patches, min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% load patches for *one profile*
%
%
% 3/7/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

tiwe_patches_paths

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)];

if merge_patches==1
    load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
else
    load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
end


%%