%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patches_all_eq08.m
%
%
%------------
% 3/29/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Make plots

clear ; close all

project_name = 'eq08'
depth_range = [60 200]

patch_size_min = 1
usetemp = 1
merge_patches = 0
min_sep = 0.15
%
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

%
h=plot_patch_locations_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep)

%%
h=plot_patch_gamma_2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)

%%
h=plot_gamma_vs_depth2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep)

%%
h=plot_gamma_vs_epsilon2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)

%%

plot_gamma_binned_eq08

%%

plot_gamma_vs_cnum_eq08

%%
