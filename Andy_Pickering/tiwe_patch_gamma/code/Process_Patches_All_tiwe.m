%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Run patch-analysis for TIWE 
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% script to run all processing !

clear ; close all

patch_size_min = 1  % min patch size
usetemp = 1

% option to use merged patches
merge_patches = 0 ;
min_sep = 0.15 ;

% range of casts to process
cnums_to_do = 1:4000
% [2836:3711] == ydays 324-327 for TIWE

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

%
% * don't need to run FindPatches again for merged **
FindPatches_cham_Raw('tiwe',patch_size_min,usetemp,cnums_to_do)
%

if merge_patches==1
    merge_patches_tiwe(patch_size_min,usetemp,...
        merge_patches,min_sep)
end

%
Compute_N2_dTdz_patches_eachcast('tiwe',patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
add_binned_to_patches_tiwe(patch_size_min,usetemp,...
    merge_patches,min_sep, cnums_to_do)
%
Run_tiwe_AP_forPatches(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)

%
add_patch_chi_eps_to_patches_each_profile('tiwe',patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)

%
combine_patch_profiles('tiwe',patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)

%
add_R2_to_patches('tiwe',patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)

%%