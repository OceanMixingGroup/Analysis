function avg = load_avg_patch_eq14(cnum,patch_size_min,usetemp,...
    merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Load chameleon profiles processed over just patches
% These are made in run_eq14_for_PATCHES.m
%
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

eq14_patches_paths

% folder for chameleon data processed over patches (run_eq14_for_Patches.m)
if merge_patches==1
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_merged_minsep_' num2str(min_sep*100)])
else
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
end

fname=['eq14_' sprintf('%04d',cnum) '.mat'];

load(fullfile(data_dir,fname))


%%