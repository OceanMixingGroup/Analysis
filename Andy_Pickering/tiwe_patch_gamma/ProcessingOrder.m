%
% order of processing scripts
%
%
%%

FindPatches_tiwe_Raw
Compute_N2_dTdz_patches_tiwe_eachcast.m
add_binned_to_patches
Run_tiwe_AP_forPatches
add_patch_chi_eps_to_patches_tiwe_each_profile
combine_patch_profiles
plot_patch_gamma_tiwe

%% script to run all processing !

clear ; close all

patch_size_min = 0.5  % min patch size
usetemp = 1

% option to use merged patches
merge_patches = 0 ;
min_sep = 0.15 ;

FindPatches_tiwe_Raw(patch_size_min,usetemp,...
    merge_patches,min_sep)

Compute_N2_dTdz_patches_tiwe_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep)

add_binned_to_patches(patch_size_min,usetemp,...
    merge_patches,min_sep)

Run_tiwe_AP_forPatches(patch_size_min,usetemp,...
    merge_patches,min_sep)

add_patch_chi_eps_to_patches_tiwe_each_profile(patch_size_min,...
    usetemp,merge_patches,min_sep)

combine_patch_profiles(patch_size_min,usetemp,...
    merge_patches,min_sep)

%%