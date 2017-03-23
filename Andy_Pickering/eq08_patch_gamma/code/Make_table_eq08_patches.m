%~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_table_eq08_patches.m
%
%
%-----------------
% 3/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all ; clc

% make table of Npatches, avg. gamma, etc for different OT size and
merge_patches=0;min_sep=0;
and=' & ';
clc

for whp=[1 2 3]
    
    switch whp
        case 1
            patch_size_min=0.4;
            usetemp=1;
        case 2
            patch_size_min=0.75;
            usetemp=1;
        case 3
            patch_size_min=1;
            usetemp=1;
    end
    
    patches = load_eq08_patches_comb(patch_size_min, usetemp, merge_patches, min_sep);
    minR2=0;
    depth_range = [60 200];
    id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2) & patches.R2>minR2) ;
    
    disp([num2str(patch_size_min) and num2str(usetemp) and num2str(minR2) and num2str(roundx(nanmedian(patches.gam_bin(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_line(id)),2)) and num2str(roundx(nanmedian(patches.gam_line_fit(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_bulk(id)),2)) and num2str(length(patches.p1(id))) ...
        ' \\' ])
    
    disp('\hline')
    minR2=0.5;
    id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2) & patches.R2>minR2) ;    
    
    disp([num2str(patch_size_min) and num2str(usetemp) and num2str(minR2) and num2str(roundx(nanmedian(patches.gam_bin(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_line(id)),2)) and num2str(roundx(nanmedian(patches.gam_line_fit(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_bulk(id)),2)) and num2str(length(patches.p1(id))) ...
        ' \\' ])    
    
    disp('\hline')
    
end

%%