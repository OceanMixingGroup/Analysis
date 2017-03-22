function [] = combine_patch_profiles_eq14(patch_size_min,...
    usetemp,merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% combine_patch_profiles_eq14.m
%
% Combine patch data from each profile into a single structure.
%
% *add field Npatches for each cnum
%
%-----------
% 2/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%clear ; close all

% patch options
%patch_size_min = 0.25  % min patch size
%usetemp = 1

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]


% set paths
eq14_patches_paths

ip=0;
hb=waitbar(0,'combine_patch_profiles_eq14');

for cnum=1:3200
    waitbar(cnum/3200,hb)
    try
        
        % load the patches for this profile
        clear patch_data patches
        if merge_patches==1
            load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
        else
            load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
        end
        
        ip=ip+1;
        
        if ip==1
            patch_all=patches;
        else
            % add each field on to patch_all
            fnames=fieldnames(patches);
            for ivar=1:length(fnames)
                patch_all.(fnames{ivar}) = [patch_all.(fnames{ivar})(:) ; patches.(fnames{ivar})(:)] ;
            end
            
        end % ip==1
        
        if length(patch_all.cnum)~=length(patch_all.gam_range)
            disp(['uhoh ' num2str(cnum)])
        end
        
    end % try
end % cnum

delete(hb)
%%

clear patches
patches=patch_all; clear patch_all
patches.MakeInfo = ['Made ' datestr(now) ' w/ combine_patch_profiles_eq14.m']

%%

% save combined structure
if merge_patches==1
    save( fullfile( analysis_dir, project_long, 'data',...
        [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
else
    save( fullfile( analysis_dir, project_long, 'data',...
        [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )
end
%%