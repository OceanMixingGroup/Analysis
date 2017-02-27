%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_binned_to_patches_eq14.m
%
% Add binned (standard 1m avg data) chameleon data to patch data for each
% profile. Binned data are interpolated to patch locations. Patch data for
% each profile are re-saved with added binned values.
%
%
%----------------
% 2/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Add 1m binned data (interpolated to patch locations)

clear ; close all

% patch parameters
patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

eq14_patches_paths
%

% load binned chameleon data (structure containing all profiles)
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

hb=waitbar(0)

for cnum=1:3200
    
    waitbar(cnum/4000,hb)
    
    clear patches Npatches
    
    try
        
        % load patch data for this profile
        clear patch_data patches
        load(fullfile(save_dir_patch,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']))
        
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
        %    patches.drhodz_bin(ip)= patches.n2_bin(ip) * (nanmean(cham.SIGMA(:,Icham))+1000) / -9.81 ;
        
        %     if log10(patches.eps_bin(ip))>-8.5
        patches.gam_bin=ComputeGamma(patches.n2_bin,patches.dtdz_bin,patches.chi_bin,patches.eps_bin);
        %     end
        
        % re-save profile
        save(fullfile(save_dir_patch,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        
    end % try
end % cnum
delete(hb)

%%