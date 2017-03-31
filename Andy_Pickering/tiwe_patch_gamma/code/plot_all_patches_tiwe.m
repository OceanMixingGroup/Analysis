%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_all_patches_tiwe.m
%
%
%------------
% 3/29/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% non-patch plots

clear ; close all

plot_TIWE_avg_comb

h=plot_gamma_binned

Plot_T_S_raw_tiwe

%% make plots

clear ; close all

for  whcase = 1
    
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
    h=plot_patch_locations_tiwe(patch_size_min,usetemp,...
        merge_patches,min_sep)
       
    %
    h=plot_patch_gamma_2X2('tiwe',patch_size_min,usetemp,...
    merge_patches,min_sep,[80 200])

    %
    h=plot_gamma_vs_yday(patch_size_min,usetemp,...
        merge_patches,min_sep,[80 200])

    %
    h=compare_patches_tiwe_AP_Bill(patch_size_min,usetemp,...
        merge_patches,min_sep,0.0)
    
    %
    h= plot_gamma_vs_depth2X2('tiwe',patch_size_min,usetemp,...
        merge_patches,min_sep)
    
    %
    h=plot_gamma_vs_epsilon2X2('tiwe',patch_size_min,usetemp,...
        merge_patches,min_sep,[80 200])
    
end
%%