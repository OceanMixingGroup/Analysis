%
%
% eq14_patches_paths.m
%
%
%
%%

clear analysis_dir project_long project_short

analysis_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering' ;
project_long      = 'eq14_patch_gamma' ;
project_short= 'eq14' ;
save_dir_patch = fullfile(analysis_dir,project_long,'data','patches') ;

save_dir_avg   = fullfile(analysis_dir,project_long,'data','avg') ;

save_dir_avg_patch = fullfile(analysis_dir,project_long,'data','avg_patch') ;

save_dir_cal='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal' ;

fig_dir = fullfile( analysis_dir,project_long,'figures');

%%

