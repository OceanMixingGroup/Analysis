function []=add_binned_to_patches_eq08(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_binned_to_patches_eq08.m
%
% Add binned (standard 1m avg data) chameleon data to patch data for each
% profile. Binned data are interpolated to patch locations. Patch data for
% each profile are re-saved with added binned values.
%
%
%
%----------------
% 3/21/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Add 1m binned data (interpolated to patch locations)


ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]

% set paths
eq08_patches_paths

% load binned chameleon data (structure containing all profiles)
%***
%load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_1mavg_combined.mat')
%***

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/


hb=waitbar(0)
ic=0;
for cnum=cnums_to_do %2836:3711 %1:4000
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    
    clear patches Npatches
    
    try
        
        % load the patches for this profile
        if merge_patches==1
        load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )    
        else
        load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
        end
        
        % make empty arrays for the binned data
        Npatches = length(patches.cnum) ;
        patches.gam_bin = nan*ones(size(patches.p1)) ;
        patches.n2_bin  = nan*ones(size(patches.p1)) ;
        patches.dtdz_bin= nan*ones(size(patches.p1)) ;
        patches.chi_bin = nan*ones(size(patches.p1)) ;
        patches.eps_bin = nan*ones(size(patches.p1)) ;
        patches.drhodz_bin=nan*ones(size(patches.p1)) ;
        
        % get index for binned profile
        clear Icham pbin pmn
        Icham= find(cham.cnum==cnum) ;
        pbin = cham.P;
        pmn  = nanmean([patches.p1 patches.p2],2) ;
        
        % interp binned data to the patch locations (mean p of each patch)
        
        clear ig
        ig=find(~isnan(pbin)) ;
        patches.n2_bin    = interp1(pbin(ig) , cham.N2(ig,Icham) , pmn);
        patches.dtdz_bin  = interp1(pbin(ig) , cham.DTDZ(ig,Icham) , pmn);
        patches.chi_bin   = interp1(pbin(ig) , cham.CHI(ig,Icham) , pmn);
        patches.eps_bin   = interp1(pbin(ig) , cham.EPSILON(ig,Icham) , pmn);
        
        clear ib
        ib = find(patches.eps_bin<0.4e-9);
        patches.eps_bin(ib) = nan ;

        patches.gam_bin=ComputeGamma(patches.n2_bin,patches.dtdz_bin,patches.chi_bin,patches.eps_bin);
        
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