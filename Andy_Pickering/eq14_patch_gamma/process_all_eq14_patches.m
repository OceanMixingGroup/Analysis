
%%

clear all ; close all

patch_size_min = 1
usetemp = 1

merge_patches = 0
min_sep = 0.15

%cnums_to_do = [4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
%        1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];
    
cnums_to_do = [1:5:3100];
%

FindPatches_EQ14_RAw(patch_size_min,usetemp, cnums_to_do)

%

if merge_patches==1
merge_patches_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)
end

%
Compute_N2_dTdz_patches_eq14_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)

%

add_binned_to_patches_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
run_eq14_for_PATCHES(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)

%
add_patch_chi_eps_to_patches_tiwe_each_profile_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
combine_patch_profiles_eq14(patch_size_min,...
    usetemp,merge_patches,min_sep)

%
add_R2_to_patches_eq14(patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)

%%


