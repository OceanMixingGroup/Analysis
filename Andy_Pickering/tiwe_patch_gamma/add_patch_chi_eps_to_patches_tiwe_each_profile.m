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

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1
save_dir_patch='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/'

% folder for chameleon data processed over patches (Run_tiwe_AP_forPatches.m)
data_dir=fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/avg_patch',['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

hb=waitbar(0,'compiling patch data from all profiles');

for cnum=1:4000
    
    waitbar(cnum/4000,hb)
    
    try
        % load the patches for this profile
        clear patches
        load(fullfile(save_dir_patch,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
        
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
            ib=find(log10(patches.eps)<-8.5);
            patches.eps(ib)=nan;
            
            % compute gamma for each patch
            patches.gam_range = ComputeGamma(patches.n2_range, patches.dtdz_range, patches.chi, patches.eps);
            patches.gam_line  = ComputeGamma(patches.n2_line , patches.dtdz_line , patches.chi, patches.eps);
            patches.gam_bulk  = ComputeGamma(patches.n2_bulk , patches.dtdz_bulk , patches.chi, patches.eps);
            patches.gam_bulk_2= ComputeGamma(patches.n2_bulk_2, patches.dtdz_bulk, patches.chi, patches.eps);
            patches.gam4 = ComputeGamma(patches.n4, patches.dtdz_line, patches.chi , patches.eps );
            
        end
        % re-save profile
        save(fullfile(save_dir_patch,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        
    end % try
    
end % cnum
delete(hb)

%%
