%~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_table_tiwe_patches.m
%
% Make latex table of Npatches, avg. gamma, etc for different parameters
%
%
%-----------------
% 3/10/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all ; clc

merge_patches=0;min_sep=0;
depth_range=[80 200];
and=' & ';
clc

for whp=[1  3]%[1 3]
    
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
        case 4
            patch_size_min=1;
            usetemp=0;
    end
    
    patches = load_tiwe_patches_comb(patch_size_min, usetemp, merge_patches, min_sep);
    minR2=0;
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

%% Same, but for only a specified yday range

clear ; close all ; clc

merge_patches=0;min_sep=0;
depth_range=[80 200];
yday_range=[324 327];
and=' & ';
clc

for whp=[1  3]%[1 3]
    
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
        case 4
            patch_size_min=1;
            usetemp=0;
    end
    
    patches = load_tiwe_patches_comb(patch_size_min, usetemp, merge_patches, min_sep);
    minR2=0;
    id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2) & patches.R2>minR2 ...
        & patches.yday>=yday_range(1) & patches.yday<=yday_range(2) ) ;
    
    disp([num2str(patch_size_min) and num2str(usetemp) and num2str(minR2) and num2str(roundx(nanmedian(patches.gam_bin(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_line(id)),2)) and num2str(roundx(nanmedian(patches.gam_line_fit(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_bulk(id)),2)) and num2str(length(patches.p1(id))) ...
        ' \\' ])
    
    disp('\hline')
    minR2=0.5;
    id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2) & patches.R2>minR2 ...
                & patches.yday>=yday_range(1) & patches.yday<=yday_range(2) )  ;
    
        disp([num2str(patch_size_min) and num2str(usetemp) and num2str(minR2) and num2str(roundx(nanmedian(patches.gam_bin(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_line(id)),2)) and num2str(roundx(nanmedian(patches.gam_line_fit(id)),2)) ...
        and num2str(roundx(nanmedian(patches.gam_bulk(id)),2)) and num2str(length(patches.p1(id))) ...
        ' \\' ])

    
    disp('\hline')
    
end
%%