function [] = add_patch_chi_eps_to_patches_tiwe_each_profile_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_patch_chi_eps_to_patches_tiwe_each_profile_eq08.m
%
% Add chi and epsilon values computed for patches (in
% run_eq14_for_PATCHES.m) to our patches structure for each profile.
%
% Modified from add_patch_chi_eps_to_patches_tiwe_each_profile.m
%
%
%-----------
% 3/22/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% set paths
eq08_patches_paths

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]


% folder for chameleon data processed over patches (run_eq14_for_Patches.m)
if merge_patches==1
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_merged_minsep_' num2str(min_sep*100)])
else
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
end

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

hb=waitbar(0,'add_patch_chi_eps');
ic=0;
for cnum=cnums_to_do%1:3200
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    
    try
        % load the patches for this profile
        clear patches
        if merge_patches==1
            load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
        else
            load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
        end
        
        % add empty arrays for chi and eps
        patches.eps=nan*ones(size(patches.p1));
        patches.chi=nan*ones(size(patches.p1));
        
        patches.gam_line = nan*ones(size(patches.p1));
        patches.gam_bulk = nan*ones(size(patches.p1));
        patches.gam4 = nan*ones(size(patches.p1));
        
        % re-save profile
        if merge_patches==1
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
        else
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        end
        
    end
    
    try
        % load the processed profile w/ chi and eps for each patch ('avg')
        clear avg
        fname=['eq08' sprintf('%04d',cnum) '_avg.mat'];
        load(fullfile(data_dir,fname))
        
        if length(patches.p1)==length(avg.CHI)
            
            patches.chi = avg.CHI(:);
            patches.eps = avg.EPSILON(:);
            
            clear ib
            ib=find(log10(patches.eps)<-8.5);
            patches.eps(ib)=nan;
            
            % compute gamma for each patch
            %patches.gam_range = ComputeGamma(patches.n2_range, patches.dtdz_range, patches.chi, patches.eps);
            patches.gam_line  = ComputeGamma(patches.n2_line , patches.dtdz_line , patches.chi, patches.eps);
            patches.gam_bulk  = ComputeGamma(patches.n2_bulk , patches.dtdz_bulk , patches.chi, patches.eps);
            %patches.gam_bulk_2= ComputeGamma(patches.n2_bulk_2, patches.dtdz_bulk, patches.chi, patches.eps);
            patches.gam4 = ComputeGamma(patches.n4, patches.dtdz_line, patches.chi , patches.eps );
            
        end
        
        % re-save profile
        if merge_patches==1
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
        else
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        end
        
    catch
        disp('error~!!!')
    end % try
    
    
end % cnum
delete(hb)

%%
