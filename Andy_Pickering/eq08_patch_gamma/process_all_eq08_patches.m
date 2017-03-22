
%%
clear ; close all

patch_size_min = 0.4
usetemp = 1
merge_patches = 0
min_sep = 0.15

cnums_to_do = [500:600] ;

%%
FindPatches_EQ08_Raw(0.4,1,cnums_to_do )

%%
%

% if merge_patches==1
% merge_patches_eq14(patch_size_min,usetemp,...
%     merge_patches,min_sep)
% end

%%
Compute_N2_dTdz_patches_eq08_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep, cnums_to_do)

%%
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