%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% merge_patches.m
%
% Try to merge patches separated < min distance, as done in Smyth et al
% 2001
%
%
%-------------
% 3/1/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load (un-merged) patches for a profile
load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/tiwe_raw_patches_minOT_15_usetemp_1_cnum_1902.mat')
%
%patch_data = patch_data(1:5,:)

Np=size(patch_data,1)

min_sep = 0.15;

%new_patches=[patch_data(1,:)]
new_patches = patch_data
%
% for ip=2:Np
%     ip
%     if ip<Np
ip=2
while ip<Np
    ip
    if ( new_patches(ip,2) - new_patches(ip-1,3)  ) < min_sep % separated by < min_sep
        disp('merging patches')
        clear merged_patch
        % merge patch w/ one before it
        new_patches(ip-1,:)=[ new_patches(ip,1) new_patches(ip-1,2) new_patches(ip,3) nan nan nan new_patches(ip,7) ];
        % remove the old patch we merged
        new_patches = [ new_patches(1:(ip-1),:) ; new_patches( (ip+1):Np,:) ];
        Np=Np-1        
    else
        ip=ip+1
    end
end
size(new_patches)
Np
new_patches;
%end

size(new_patches)

%%

figure(1);clf
plot(1,patch_data(:,2),'kp')
hold on
plot(1,patch_data(:,3),'rd')
axis ij
plot(2,new_patches(:,2),'kp')
hold on
plot(2,new_patches(:,3),'rd')
xlim([0 3])

%%

[patch_data(1:10,2) new_patches(1:10,2)]

