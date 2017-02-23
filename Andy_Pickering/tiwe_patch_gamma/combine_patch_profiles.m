%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% combine_patch_profiles.m
%
%
% 2/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1
save_dir_patch='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/'

ip=0;
hb=waitbar(0)
for cnum=1:4000
    waitbar(cnum/4000,hb)
    try
        
        % load the patches for this profile
        clear patches
        load(fullfile(save_dir_patch,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
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
        
    end % try
end % cnum

delete(hb)

clear patches
patches=patch_all; clear patch_all

% save combined structure
save( fullfile( '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/',...
    ['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%%