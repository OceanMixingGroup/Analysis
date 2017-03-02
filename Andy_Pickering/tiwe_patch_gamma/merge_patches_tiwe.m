%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% merge_patches_tiwe.m
%
% Merge patches that are separated by less than a minimum distance, as done
% in Smyth et al 2001.
%
%------------
% 3/1/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

% patch options
patch_size_min = 0.15  % min patch size
usetemp = 1

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]

% min separation for merging
min_sep = 0.15;

% set paths
tiwe_patches_paths

hb=waitbar(0)
warning off

for cnum=1:4000
    waitbar(cnum/4000,hb)
    try
        
        % load patch data for this profile
        clear patch_data 
        load(fullfile(save_dir_patch,ot_dir,[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']))
        
        clear new_patches ip Np
        new_patches = patch_data;
        Np=size(new_patches,1);
        %
        ip=2;
        while ip<Np
            ip;
            if ( new_patches(ip,2) - new_patches(ip-1,3)  ) < min_sep % separated by < min_sep
                %disp('merging patches')
                clear merged_patch
                % merge patch w/ one before it
                new_patches(ip-1,:)=[ new_patches(ip,1) new_patches(ip-1,2) new_patches(ip,3) nan nan nan new_patches(ip,7) ];
                % remove the old patch we merged
                new_patches = [ new_patches(1:(ip-1),:) ; new_patches( (ip+1):Np,:) ];
                Np=Np-1;
            else
                ip=ip+1;
            end
        end 
        
        clear patch_data
        patch_data = new_patches ;
        
        % save profile
        save(fullfile(save_dir_patch,ot_dir,[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']), 'patch_data' )
    end % try
end % cnum

delete(hb)
warning on
%%