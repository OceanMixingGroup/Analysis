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

for whcase = 9:11
    
    switch whcase
        
        case 1
            
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 0 ; Pmin  = 0 ;
            
        case 2
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 3
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 10 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 4
            
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 10 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 50 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 5
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 0 ; Pmin  = 0 ;
            
        case 6
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 7
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 10 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 8
            this_proj = 'eq08'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 10  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 50 ;
            cnums_to_get = 200:2700;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
            
            
            %~~~~~~~~~~~ EQ14 ~~~~~~~~~~~~~~~~~
            
        case 9
            
            this_proj = 'eq14'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 7  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = get_cham_cnums_eq14 ;
            screen_chi= 1 ; screen_ml = 0 ; Pmin  = 0 ;
            
        case 10
            this_proj = 'eq14'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 7  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 2 ;
            cnums_to_get = get_cham_cnums_eq14 ;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 11
            this_proj = 'eq14'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 7  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 10 ;
            cnums_to_get = get_cham_cnums_eq14;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
        case 12
            
            this_proj = 'eq14'
            eval([this_proj '_patches_paths'])
            
            Params.gamma    = 0.2;
            Params.fmax     = 7  ;
            Params.z_smooth = 1 ;
            Params.resp_corr= 0  ;
            Params.fc       = 99 ;
            
            dz = 50 ;
            cnums_to_get = get_cham_cnums_eq14;
            screen_chi= 1 ; screen_ml = 1 ; Pmin  = 20 ;
            
            
    end % switch
    
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml) ;
    
    addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/
    
    save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
    
    save( fullfile(analysis_dir,project_short,'data',save_name), 'chipod','cham')
    
end % whcase
