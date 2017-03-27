function [] = add_binned_to_patches(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_binned_to_patches.m
%
% Add binned (standard 1m avg data) chameleon data to patch data for each
% profile. Binned data are interpolated to patch locations. Patch data for
% each profile are re-saved with added binned values.
%
%
%----------------
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Add 1m binned data (interpolated to patch locations)

% set paths
eval([project_name '_patches_paths'])

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]

% load binned chameleon data (structure containing all profiles)
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

%addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

hb=waitbar(0,'add binned to patches')

for cnum=1:4000
    
    waitbar(cnum/4000,hb)
    
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
        Icham= find(cham.castnumber==cnum) ;
        pbin = cham.P(:,Icham);
        pmn  = nanmean([patches.p1 patches.p2],2) ;
        
        % interp binned data to the patch locations (mean p of each patch)
        
        clear ig
        ig=find(~isnan(pbin)) ;
        patches.n2_bin    = interp1(pbin(ig) , cham.N2(ig,Icham) , pmn);
        patches.dtdz_bin  = interp1(pbin(ig) , cham.DTDZ(ig,Icham) , pmn);
        patches.chi_bin   = interp1(pbin(ig) , cham.CHI(ig,Icham) , pmn);
        patches.eps_bin   = interp1(pbin(ig) , cham.EPSILON(ig,Icham) , pmn);
        %patches.drhodz_bin(ip)= patches.n2_bin(ip) * (nanmean(cham.SIGMA(:,Icham))+1000) / -9.81 ;
        
        clear ib
        ib = find( log10(patches.eps_bin)<(-8.5) );
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