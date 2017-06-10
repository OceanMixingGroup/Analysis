%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% make_combined_data_files.m
%
% Getting and averaging/combining all the chipod and cham profiles takes a
% long time, so here i'll do that and save the results for loading .
%
%
%~~~~~~~~~~~~~~~~
% 6/9/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all


addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

for whcase=5:8
    
    switch whcase
        
        case 1
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 0 ; Pmin  = 0 ;
            
        case 2
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 3
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 10 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 4
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 50 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 5
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 0 ; Pmin  = 0 ;
            
        case 6
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 7
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 10 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 8
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            eq08_patches_paths
            dz = 50 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
    end % switch
    
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
    
    addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/
    
    save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
    
    save( fullfile(analysis_dir,project_short,'Data',save_name), 'chipod','cham')
    
end % whcase
