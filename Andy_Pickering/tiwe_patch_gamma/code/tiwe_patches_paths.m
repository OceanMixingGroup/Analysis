%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% tiwe_patches_paths.m
%
%
%
%-----------------
% 4/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

analysis_dir  = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering' ;
project_long  = 'tiwe_patch_gamma' ;
project_short = 'tiwe' ;

save_dir_patch = fullfile(analysis_dir,project_long,'data','patches') ;

save_dir_avg   = fullfile(analysis_dir,project_long,'data','avg') ;

save_dir_avg_patch = fullfile(analysis_dir,project_long,'data','avg_patch') ;

save_dir_cal = fullfile(analysis_dir,project_long,'data','cham_proc','cal') ;

fig_dir = fullfile(analysis_dir,project_long,'figures') ;

% chameleon 1m avg processed
path_cham_avg       = fullfile(analysis_dir,project_long,'data','cham_proc','avg');

% chipod method applied to 1m binned  data
path_chipod_bin     = fullfile(analysis_dir,project_long,'data','cham_proc','chipod_method_bin');

%%

