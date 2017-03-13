
%%

clear all ; close all

patch_size_min = 0.4
usetemp=1

merge_patches=1
min_sep=0.15

%%

FindPatches_EQ14_RAw(patch_size_min,usetemp)

%%

merge_patches_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)

%%

Compute_N2_dTdz_patches_eq14_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep)

%

add_binned_to_patches_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)
%

run_eq14_for_PATCHES(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
add_patch_chi_eps_to_patches_tiwe_each_profile_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)

%

combine_patch_profiles_eq14(patch_size_min,...
    usetemp,merge_patches,min_sep)

%%



