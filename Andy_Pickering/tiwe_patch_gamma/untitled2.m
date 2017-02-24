%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1
save_dir_patch='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/'

ip=0;
hb=waitbar(0)

bad=[]

for cnum=1:3400
    waitbar(cnum/4000,hb)
    try
        
        % load the patches for this profile
        clear patches
        load(fullfile(save_dir_patch,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
%        ip=ip+1;
      
if length(patches.cnum)~=length(patches.gam_range)
    bad=[bad cnum];
end

    end
end
        %%