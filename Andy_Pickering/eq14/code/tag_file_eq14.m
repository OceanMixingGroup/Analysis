% tag_file_eq14.m
% tag flags data that should be set to NaN
% integer part of the number means number of cast and 
% decimal part is the depth in meters, e.g.
% 138.136 138.137 means cast 136 depth range 136:137m
% should be set to NaNs; .999 means to the bottom
% t - temperature; c - conductivity; e - epsilon

% clear
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot

% load('~/GDrive/data/eq14/chameleon/processed/attempting_to_fix_epsilon/sum/EQ14_sum.mat')


% 
% figure(1)
% clf
% 
% ax(1) = subplot(311);
% pcolor(cham.castnumber,cham.depth,cham.SAL)
% shading flat
% axis ij
% % colormap(ax(1),flipud(lbmap(32,'RedBlue')))
% colormap(ax(1), cbrewer('seq','YlGnBu',32))
% colorbar
% caxis([34.2 35.2])
% ylabel({'SAL [psu]';'depth [m]'})
% 
% % ax(2) = subplot(312);
% ax(1) = subplot(311);
% pcolor(cham.castnumber,cham.depth,cham.T1)
% shading flat
% axis ij
% colormap(ax(1), flipud( cbrewer('div','RdYlBu',28)))
% colorbar
% caxis([12 26])
% ylabel({'T1 [^oC]';'depth [m]'})
% 
% ax(3) = subplot(313);
% pcolor(cham.castnumber,cham.depth,cham.T2)
% shading flat
% axis ij
% colormap(ax(3), flipud(cbrewer('div','RdYlBu',28)))
% colorbar
% caxis([12 26])
% ylabel({'T2 [^oC]';'depth [m]'})
% 
% ax(3) = subplot(313);
% pcolor(cham.castnumber,cham.depth,cham.COND)
% shading flat
% axis ij
% colormap(ax(3), flipud(cbrewer('div','RdYlBu',28)))
% colorbar
% caxis([4 6])
% ylabel({'COND [seimans m^{-1}]';'depth [m]'})
% 
% ax(2) = subplot(312);
% pcolor(cham.castnumber,cham.depth,log10(cham.EPSILON1))
% shading flat
% axis ij
% colormap(ax(2), jet)
% colorbar
% caxis([-9 -4])
% ylabel({'EPSILON1 [W kg^{-1}]';'depth [m]'})
% 
% ax(3) = subplot(313);
% pcolor(cham.castnumber,cham.depth,log10(cham.EPSILON2))
% shading flat
% axis ij
% colormap(ax(3), jet)
% colorbar
% caxis([-9 -4])
% ylabel({'EPSILON2 [W kg^{-1}]';'depth [m]'})
% 
% linkaxes(ax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
tag.SAL = ...
  [0004.000 0004.999;...
   0058.056 0058.057;...
   0416.049 0416.049;...
   0454.000 0454.999;...
   0455.034 0455.036;...
   0607.042 0607.044;...
   0649.018 0649.018;...
   0767.086 0767.999;...
   0785.123 0785.123;...
   0800.000 0800.012;...
   0801.000 0801.120;...
   0802.010 0802.012;...
   0834.016 0834.018;...
   0884.015 0884.015;...
   0897.022 0897.022;...
   0938.007 0938.014;...
   0939.016 0939.016;...
   1005.000 1005.999;...
   1033.000 1033.030;...
   1070.000 1074.999;...
   1181.099 1181.100;...
   1257.000 1257.999;...
   1270.000 1270.999;...
   1304.022 1304.025;...
   1340.022 1340.024;...
   1515.013 1515.015;...
   1578.000 1578.999;...
   1642.000 1642.999;...
   1675.052 1675.056;...
   1745.000 1745.011;...
   1832.016 1832.021;...
   1856.121 1856.136;...
   1895.069 1895.074;...
   1941.039 1941.044;...
   1947.032 1947.037;...
   2004.000 2004.999;...
   2072.000 2073.999;...
   2107.095 2107.999;...
   2115.011 2115.013;...
   2139.032 2139.036;...
   2282.038 2283.999;...
%    2296.056 2299.999;...
   2348.000 2348.999;...
   2371.000 2371.999;...
   2391.000 2391.999;...
   2582.000 2582.050;...
   2678.000 2678.999;...
   2695.065 2695.071;...
   2715.000 2715.999;...
   2734.019 2734.025;...
   2762.000 2762.999;...
   2800.000 2800.030;...
   2828.167 2828.172;...
   2854.000 2854.080;...
   2855.113 2855.119;...
   2935.140 2935.999;...
   2953.000 2953.999;...
   2970.000 2970.999;...
   3030.000 3030.067;...
   3031.000 3031.066];
   
   

tag.T1 = ...
    [1005.000 1005.999;...
     1070.000 1074.999;...
     2282.038 2283.999;...
     2348.000 2348.999;...
     2391.000 2391.999;...
     2678.000 2678.999;...
     2762.000 2762.999;...
     2953.000 2953.999;...
     2970.000 2970.999];

tag.T2 = ...
    [1005.000 1005.999;...   
     1070.000 1074.999;...
     1642.000 1642.999;...
     2107.095 2107.999;...
     2348.000 2348.999;...
     2678.000 2678.999;...
     2970.000 2970.999];
 
 
% In cali_eq14, set the times when THETA should be based on T2 rather than
% T1 and make sure the flags are correct here
tag.THETA = ...
    [1005.000 1005.999;...
     1070.000 1074.999;...
     2348.000 2348.999;...
     2678.000 2678.999;...
     2970.000 2970.999];

   
 
tag.EPSILON1 = ...
    [0754.000 0754.999;...
     0925.153 0925.165;...
     1070.000 1074.999;...
     1240.048 1240.064;...
     1620.060 1620.999;...
     1642.000 1647.999;...
     1690.000 1690.999;...
     2107.097 2107.999;...
     2118.100 2118.999;...
     2197.000 2197.037;...
     2282.038 2283.999;...
     2285.000 2285.053;...
     2314.021 2314.107;...
     2419.030 2419.999;...
     2516.000 2517.999;...
     2556.000 2556.999;...
     2676.000 2676.999;...
     2725.000 2725.999;...
     2778.000 2778.000;...
     2813.025 2813.058;...
     2925.000 2925.999;...
     3020.000 3020.999;...
     3085.000 3085.999];
 
tag.EPSILON2 = ...
    [0011.000 0011.999;...
     0015.000 0015.999;...
     0812.072 0812.999;...
     0893.074 0893.999;...
     0925.153 0925.165;...
     1070.000 1074.999;...
     1091.000 1091.999;...
     1240.048 1240.064;...
     1642.000 1647.999;...
     2107.097 2107.999;...
     2197.000 2197.037;...
     2282.038 2283.999;...
     2285.000 2285.053;...
     2578.000 2678.999;...
     2683.000 2683.999;...
     2729.000 2729.999;...
     2735.000 2735.999;...
     2813.025 2813.058;...
     2840.000 2840.999;...
     2895.059 2895.999;...
     2925.000 2925.999;...
     2990.057 2990.089;...
     3006.000 3006.999;...
     3034.000 3039.999];
     
     
% since EPSILON is calculated as the average of both EPSILON1 and EPSILON2,
% the tag for EPSILON needs to cover both probes. CHECK THIS. There may be
% times when the better one is used in average_data_filt_gen1 instead of
% the average in which case the data need not be NaNed out.

tag.EPSILON = [tag.EPSILON1; tag.EPSILON2];