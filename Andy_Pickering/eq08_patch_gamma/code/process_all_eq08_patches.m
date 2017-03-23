
%%
clear all ; close all ; clc

patch_size_min = 1
usetemp = 1
merge_patches = 0
min_sep = 0.15

cnums_to_do = [1:5:2700] ;
%cnums_to_do = [1:5:615] ;

%%

%
FindPatches_EQ08_Raw(patch_size_min,usetemp,cnums_to_do )

%
%

% if merge_patches==1
% merge_patches_eq14(patch_size_min,usetemp,...
%     merge_patches,min_sep)
% end

%
Compute_N2_dTdz_patches_eq08_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep, cnums_to_do)
%
add_binned_to_patches_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
run_eq08_for_PATCHES(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
add_patch_chi_eps_to_patches_tiwe_each_profile_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
combine_patch_profiles_eq08(patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)
%%
add_R2_to_patches_eq08(patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)

%% Make plots

project_name = 'eq08'
depth_range = [60 200]

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

%
h=plot_patch_locations_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_patch_gamma_2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)

%
h=plot_gamma_vs_depth2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_gamma_vs_epsilon2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep)
%%