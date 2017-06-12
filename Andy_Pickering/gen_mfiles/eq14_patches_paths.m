%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% eq14_patches_paths.m
%
%
% 4/4/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear analysis_dir project_long project_short

analysis_dir   = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering' ;
project_patch   = 'eq14_patch_gamma' ;
project_short  = 'eq14' ;

save_dir_patch = fullfile(analysis_dir,project_patch,'data','patches') ;

% chameleon profiles processed for patches only
save_dir_avg_patch = fullfile(analysis_dir,project_patch,'data','avg_patch') ;

%save_dir_avg   = fullfile(analysis_dir,project_long,'data','avg') ;

% chameleon raw T,S, etc.
save_dir_cal = fullfile( analysis_dir,'eq14','data','Cham_proc_AP','cal') ;

% chipod method applied to 1m binned  data
path_chipod_bin     = fullfile( analysis_dir,'eq14','data','Cham_proc_AP','chipod_method_bin');

% chameleon 1m avg processed
path_cham_avg       = '/Users/Andy/Cruises_Research/Chipod/Cham_Eq14_Compare/data/chameleon/processed_AP_7hz/mat/';

fig_dir = fullfile( analysis_dir,project_short,'figures');

%%

