%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% combine_patch_profiles.m
%
% Combine patch data from each profile into a single structure.
%
% *add field Npatches for each cnum
%
%-----------
% 2/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1

% set paths
tiwe_patches_paths

ip=0;
hb=waitbar(0);

for cnum=1:4000
    waitbar(cnum/4000,hb)
    try
        
        % load the patches for this profile
        clear patches
        load(fullfile(save_dir_patch,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
        ip=ip+1;
        
        if ip==1
            patch_all=patches;
        else
            % add each field on to patch_all
            fnames_all=fieldnames(patch_all);
            fnames=fieldnames(patches);
            
            if length(fnames_all)~=length(fnames)
                disp('hey')
            end
            
            for ivar=1:length(fnames)
                patch_all.(fnames{ivar}) = [patch_all.(fnames{ivar})(:) ; patches.(fnames{ivar})(:)] ;
            end
            
        end % ip==1
        
        if length(patch_all.cnum)~=length(patch_all.gam_range)
            disp(['uhoh ' num2str(cnum)])
        end
        
    catch
        %disp('error')
    end % try
end % cnum

delete(hb)
%%

clear patches
patches=patch_all; clear patch_all
patches.MakeInfo = ['Made ' datestr(now) ' w/ combine_patch_profiles.m']
%%
% save combined structure
save( fullfile( analysis_dir, project, 'data',...
    [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%%