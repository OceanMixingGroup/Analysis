function avg = load_avg_patch(project_name,cnum,patch_size_min,usetemp,...
    merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%* general function
% Load chameleon profiles processed over just patches
% These are made in run_*project_name*_for_PATCHES.m
%
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

eval([project_name '_patches_paths'])

% folder for chameleon data processed over patches (run_eq14_for_Patches.m)
if merge_patches==1
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_merged_minsep_' num2str(min_sep*100)]);
else
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]);
end

fname=[project_short '_' sprintf('%04d',cnum) '.mat'];

load(fullfile(data_dir,fname))


%%