%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_patch_chi_eps_to_patches_tiwe_each_profile_eq14.m
%
% Add chi and epsilon values computed for patches (in
% run_eq14_for_PATCHES.m) to our patches structure for each profile.
%
% Modified from add_patch_chi_eps_to_patches_tiwe_each_profile.m
%
%
%-----------
% 2/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1

eq14_patches_paths
%%
%save_dir_patch='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/'

% folder for chameleon data processed over patches (run_eq14_for_Patches.m)
%data_dir=fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/avg_patch',['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
data_dir = save_dir_avg_patch

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

hb=waitbar(0,'compiling patch data from all profiles');

for cnum=1:3200
    
    waitbar(cnum/3200,hb)
    
    try
        % load the patches for this profile
        % load patch data for this profile
        clear patch_data patches
        load(fullfile(save_dir_patch,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']))

        
        % add empty arrays for chi and eps
        patches.eps=nan*ones(size(patches.p1));
        patches.chi=nan*ones(size(patches.p1));
        
        % load the processed profile w/ chi and eps for each patch ('avg')
        clear avg
        fname=['eq14_' sprintf('%04d',cnum) '.mat'];
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
        save(fullfile(save_dir_patch,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        
    end % try
    
end % cnum
delete(hb)

%%
