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

patch_size_min = .4  % min patch size
usetemp = 0

% option to use merged patches
merge_patches = 0 ;
min_sep = 0.15 ;

%%

% range of casts to process (*currently hardwired in each function)
cnum_range = [2836:3711] % ydays 324-327 for TIWE

% * don't need to run FindPatches again for merged **
FindPatches_tiwe_Raw(patch_size_min,usetemp,...
    merge_patches,min_sep)

%%

if merge_patches==1
    merge_patches_tiwe(patch_size_min,usetemp,...
    merge_patches,min_sep)
end

Compute_N2_dTdz_patches_tiwe_eachcast(patch_size_min,usetemp,...
    merge_patches,min_sep)

add_binned_to_patches(patch_size_min,usetemp,...
    merge_patches,min_sep)

Run_tiwe_AP_forPatches(patch_size_min,usetemp,...
    merge_patches,min_sep)

add_patch_chi_eps_to_patches_tiwe_each_profile(patch_size_min,...
    usetemp,merge_patches,min_sep)
%%
combine_patch_profiles(patch_size_min,usetemp,...
    merge_patches,min_sep)

%% make plots

clear ; close all

patch_size_min = 1  % min patch size
usetemp = 0

% option to use merged patches
merge_patches = 0 ;
min_sep = 0.15 ;

%%

%
h=plot_patch_locations_tiwe(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_patch_gamma_tiwe(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_gamma_vs_yday(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=compare_patches_tiwe_AP_Biill(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_gamma_vs_epsilon(patch_size_min,usetemp,...
    merge_patches,min_sep)

%
h=plot_gamma_vs_depth(patch_size_min,usetemp,...
    merge_patches,min_sep)


%% Plot gamma vs epsilon

figure(1);clf
%plot( log10(patches.gam_line),log10(patches.eps),'.')
histogram2( log10(patches.gam_line),log10(patches.eps),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}\gamma','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

%% Plot # patches vs depth?

figure(1);clf
histogram(patches.p1)
view(90, 90)
xlabel(' depth','fontsize',16)
ylabel('N patches','fontsize',16)
grid on
title(['TIWE patches - minOT=' num2str(100*patch_size_min) 'cm' ])

%%