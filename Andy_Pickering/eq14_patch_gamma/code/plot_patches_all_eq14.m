%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patches_all_eq14.m
%
%
% 3/29/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Make plots

clear ; close all

patch_size_min = 0.4
usetemp = 1
merge_patches = 0
min_sep = 0.15

eq14_patches_paths

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

%
h=plot_patch_locations_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
depth_range = [ 80 200]
h=plot_patch_gamma_2X2(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)

%
h=plot_gamma_vs_depth2X2(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_gamma_vs_epsilon2X2(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep,[0 200])

%
h=plot_gamma_vs_epsilon2X2(project_short,patch_size_min,usetemp,...
    merge_patches,min_sep,[80 200])

%%

plot_gamma_vs_cnum_eq14

%%
Plot_Processed_EQ14_cham

%%
plot_gamma_binned

%%