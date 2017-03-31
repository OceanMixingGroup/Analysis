%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% process_all_eq08_patches.m
%
% Run all patch-analysis processing for eq08
%
%
%-------------
% 3/29/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear all ; close all ; clc

patch_size_min = 1
usetemp = 1
merge_patches = 0
min_sep = 0.15

cnums_to_do = [1:5:2700] ;
%cnums_to_do = [1:5:615] ;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths
%%
FindPatches_cham_Raw(project_short,patch_size_min,usetemp,cnums_to_do )

% if merge_patches==1
% merge_patches_eq14(patch_size_min,usetemp,...
%     merge_patches,min_sep)
% end

%
Compute_N2_dTdz_patches_eachcast(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep, cnums_to_do)
%%
add_binned_to_patches_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
run_eq08_for_PATCHES(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
add_patch_chi_eps_to_patches_each_profile(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
combine_patch_profiles(project_short,patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)
%
add_R2_to_patches(project_short,patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)

%%