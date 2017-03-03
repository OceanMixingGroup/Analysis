function []=add_patch_chi_eps_to_patches_tiwe_each_profile(patch_size_min,...
    usetemp,merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_patch_chi_eps_to_patches_tiwe_each_profile.m
%
% Modified from add_patch_chi_eps_to_patches_tiwe.m
%
% Add chi and epsilon values computed for patches (in
% Run_tiwe_AP_forPatches.m) to our patches structure for each profile.
%
% Used to be part of Compute_N2_dTdz_patches_tiwe.m
%
%-----------
% 2/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%clear ; close all

% patch options
%patch_size_min = 0.4  % min patch size
%usetemp = 1

% option to use merged patches
%merge_patches = 0 ;
%min_sep = 0.15 ;

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]

% set paths
tiwe_patches_paths

% folder for chameleon data processed over patches (Run_tiwe_AP_forPatches.m)
if merge_patches==1
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_merged_minsep_' num2str(min_sep*100)])
else
    data_dir = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
end

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

hb=waitbar(0,'add_patch_chi_eps_to_patches');

% casts 2836:3711 are ydays 324-327

for cnum=2836:3711 %1:4000
    
    waitbar(cnum/4000,hb)
    
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
        
        % load the processed profile w/ chi and eps for each patch ('avg')
        clear avg
        fname=['tw91' sprintf('%04d',cnum) '_avg.mat'];
        load(fullfile(data_dir,fname))
        
        if length(patches.p1)==length(avg.CHI)
            patches.chi = avg.CHI(:);
            patches.eps = avg.EPSILON(:);
            
            clear ib
            %ib=find(log10(patches.eps)<-8.5);
            ib = find(patches.eps<0.4e-9);
            patches.eps(ib)=nan;
            
            % compute gamma for each patch
            patches.gam_range = ComputeGamma(patches.n2_range, patches.dtdz_range, patches.chi, patches.eps);
            patches.gam_line  = ComputeGamma(patches.n2_line , patches.dtdz_line , patches.chi, patches.eps);
            patches.gam_bulk  = ComputeGamma(patches.n2_bulk , patches.dtdz_bulk , patches.chi, patches.eps);
            patches.gam_bulk_2= ComputeGamma(patches.n2_bulk_2, patches.dtdz_bulk, patches.chi, patches.eps);
            patches.gam4 = ComputeGamma(patches.n4, patches.dtdz_line, patches.chi , patches.eps );
            
        end
        
        % re-save profile
        if merge_patches==1
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
        else
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        end
        
    end % try
    
end % cnum
delete(hb)

%%
