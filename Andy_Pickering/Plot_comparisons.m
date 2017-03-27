%%

project_name = 'tiwe'
%project_name = 'eq08'
%project_name = 'eq14'
patch_size_min = 0.4
depth_range = [60 200]

h=plot_patch_N2_2X2(project_name,patch_size_min,1,0,0,depth_range)

%%

h=plot_patch_Tz_2X2(project_name,patch_size_min,1,0,0,depth_range)

%%

clear

patches=load_patches_comb('tiwe',0.4,1,0,0)


