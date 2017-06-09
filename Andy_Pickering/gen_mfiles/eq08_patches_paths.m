%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% eq08_patches_paths.m
%
% Set paths for processing/analysis related to EQ08 chameleon data
%
%--------------
% 2/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

project_long      = 'eq08_patch_gamma' ;
project_short     = 'eq08' ;

analysis_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering' ;

save_dir_patch = fullfile(analysis_dir,project_long,'data','patches') ;
ChkMkDir(save_dir_patch)

save_dir_avg   = fullfile(analysis_dir,project_short,'data','avg') ;
%ChkMkDir(save_dir_avg)

save_dir_avg_patch = fullfile(analysis_dir,project_short,'data','avg_patch') ;
ChkMkDir(save_dir_avg_patch)

save_dir_cal = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc/cal/' ;

fig_dir = fullfile(analysis_dir,project_short,'figures');
ChkMkDir(fig_dir)

% chameleon 1m avg processed
%path_cham_avg       = '/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/avg/';
path_cham_avg       = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc_AP_10hz/avg/';

% chipod method applied to 1m binned  data
path_chipod_bin     = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc/chipod_method_bin/';


%%

