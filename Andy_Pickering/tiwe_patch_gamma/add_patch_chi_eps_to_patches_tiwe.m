%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% add_patch_chi_eps_to_patches_tiwe.m
%
% Add chi and epsilon values computed for patches (in
% Run_tiwe_AP_forPatches.m) to our patches structure.
%
% Used to be part of Compute_N2_dTdz_patches_tiwe.m
%
%-----------
% 2/20/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 

clear ; close all

% patch options
patch_size_min = 0.15  % min patch size
usetemp = 1

% load our 'patches' structure
load(fullfile( '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/',...
    ['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) )

%%

% folder for chameleon data processed over patches
% (Run_tiwe_AP_forPatches.m)
data_dir=fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/avg_patch',['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])

chi_all = [] ;
eps_all = [] ;
P_all = [] ;
cnum_all = [] ;

patches.eps=nan*ones(size(patches.p1));
patches.chi=nan*ones(size(patches.p1));

cnums=unique(patches.cnum);

hb=waitbar(0,'compiling patch data from all profiles');

for ic=1:length(cnums)
    waitbar(ic/length(cnums),hb)
    
    clear avg cnum
    
    cnum=cnums(ic);
    
    try
        % load the processed profile w/ chi and eps
        fname=['tw91' sprintf('%04d',cnum) '_avg.mat'];
        load(fullfile(data_dir,fname))
        
        % now loop over all patches and get average Eps and chi from those
        % depths
        clear ig
        ig=find(patches.cnum==cnum);
        for ip=1:length(ig)
            
            clear iz
            iz=isin( avg.P,[patches.p1(ig(ip)) patches.p2(ig(ip))]);
            
            if size(iz)==1 % should be 1 value for each patch
                patches.eps(ig(ip)) = nanmean(avg.EPSILON(iz)) ;
                patches.chi(ig(ip)) = nanmean(avg.CHI(iz))     ;
            end
            
        end %ip
        
    end % try
    
end % cnum
delete(hb)

%% Exclude values where epsilon is below noise floor
% note this makes significant difference in gamma median
ib=find(log10(patches.eps)<-8.5);
patches.eps(ib)=nan;

%%
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/
% use 'bulk gradient' for both
gam_range = ComputeGamma(patches.n2_range, patches.dtdz_range, patches.chi, patches.eps);
% use polyfit for both
gam_line = ComputeGamma(patches.n2_line, patches.dtdz_line, patches.chi, patches.eps);
% use bulk for both
gam_bulk = ComputeGamma(patches.n2_bulk, patches.dtdz_bulk, patches.chi, patches.eps);
%
gam_bulk_2 = ComputeGamma(patches.n2_bulk_2, patches.dtdz_bulk, patches.chi, patches.eps);
%
gam4 = ComputeGamma(patches.n4, patches.dtdz_line, patches.chi , patches.eps );

% 
patches.gam_range = gam_range;
patches.gam_line  = gam_line;
patches.gam_bulk  = gam_bulk;
patches.gam_bulk_2 = gam_bulk_2;
patches.gam4 = gam4;
%%
save( fullfile( '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/',...
    ['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%%