%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% order of processing scripts:
%
% FindPatches_tiwe_Raw
% Compute_N2_dTdz_patches_tiwe_eachcast.m
% add_binned_to_patches
% Run_tiwe_AP_forPatches
% add_patch_chi_eps_to_patches_tiwe_each_profile
% combine_patch_profiles
% plot_patch_gamma_tiwe
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
cnums_to_do = [2836:3711] % ydays 324-327 for TIWE

% * don't need to run FindPatches again for merged **
FindPatches_cham_Raw('tiwe',patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%

if merge_patches==1
    merge_patches_tiwe(patch_size_min,usetemp,...
        merge_patches,min_sep)
end

%
Compute_N2_dTdz_patches_tiwe_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%
add_binned_to_patches(patch_size_min,usetemp,...
    merge_patches,min_sep, cnums_to_do)

Run_tiwe_AP_forPatches(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
add_patch_chi_eps_to_patches_tiwe_each_profile(patch_size_min,...
    usetemp,merge_patches,min_sep)

%
combine_patch_profiles(patch_size_min,usetemp,...
    merge_patches,min_sep)

add_R2_to_patches_tiwe(patch_size_min,...
    usetemp,merge_patches,min_sep)

%% non-patch plots

plot_TIWE_avg_comb

h=plot_gamma_binned

Plot_T_S_raw_tiwe

%% make plots

clear ; close all

for  whcase = 1:3
    
    close all
    
    switch whcase
        case 1
            patch_size_min = 0.4 ; usetemp = 1 ;merge_patches = 0 ; min_sep = 0.15 ;
        case 2
            patch_size_min = 0.75 ; usetemp = 1 ;merge_patches = 0 ; min_sep = 0.15 ;
        case 3
            patch_size_min = 1 ; usetemp = 1 ;merge_patches = 0 ; min_sep = 0.15 ;
    end
    
    %
    
    %
    h=plot_patch_locations_tiwe(patch_size_min,usetemp,...
        merge_patches,min_sep)
    
    %
    h=plot_patch_gamma_tiwe(patch_size_min,usetemp,...
        merge_patches,min_sep)
    
    %
    h=plot_patch_gamma_tiwe_2X2(patch_size_min,usetemp,...
        merge_patches,min_sep)
    %
    %
    h=plot_gamma_vs_yday(patch_size_min,usetemp,...
        merge_patches,min_sep)
    
    %
    h=compare_patches_tiwe_AP_Bill(patch_size_min,usetemp,...
        merge_patches,min_sep,0.0)
    
    %
    h=plot_gamma_vs_epsilon2X2(patch_size_min,usetemp,...
        merge_patches,min_sep)
    
    %
    h=plot_gamma_vs_depth2X2(patch_size_min,usetemp,...
        merge_patches,min_sep)
    
end
%%