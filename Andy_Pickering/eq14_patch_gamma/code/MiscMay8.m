%
% Plot N2 vs dTdz, colored by epsilon?
%
%
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;
dz = 2 ;
cnums_to_get = get_cham_cnums_eq14;
%cnums_to_get = 1000:2000;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

screen_chi=1
Pmin=20
screen_ml=0

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);


%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
scatter3(log10(chipod.N2(:)), log10(chipod.Tz(:)), log10(chipod.eps(:)),30,log10(chipod.eps(:)),'filled')
%scatter3(log10(cham.N2(:)), log10(cham.Tz(:)), log10(cham.eps(:)),30,log10(cham.eps(:)),'filled')
xlabel('log_{10}[N^2]','fontsize',16)
ylabel('log_{10}[Tz]','fontsize',16)
zlabel('log_{10}[\epsilon]','fontsize',16)
zlim([-8.5 -4])
cb=colorbar;
cb.Label.String = 'log_{10}[\epsilon]';
cb.FontSize=16;
caxis([-9 -4])
title('eq14 \chi pod')
view(0,90)

subplot(212)
%scatter3(log10(chipod.N2(:)), log10(chipod.Tz(:)), log10(chipod.eps(:)),30,log10(chipod.eps(:)),'filled')
scatter3(log10(cham.N2(:)), log10(cham.Tz(:)), log10(cham.eps(:)),30,log10(cham.eps(:)),'filled')
xlabel('log_{10}[N^2]','fontsize',16)
ylabel('log_{10}[Tz]','fontsize',16)
zlabel('log_{10}[\epsilon]','fontsize',16)
zlim([-8.5 -4])
cb=colorbar;
cb.Label.String = 'log_{10}[\epsilon]';
cb.FontSize=16;
caxis([-9 -4])
title('eq14 Chameleon')
view(0,90)


%%
print( fullfile( fig_dir, 'eq14_cham_chi_N2_vs_Tz_Epscol'),'-dpng')

%%

%%

figure(1);clf
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
caxis([-9 -4])

[ir,ic]=find(log10(chipod.eps)>-4);
hold on
plot(chipod.cnum(ic),chipod.P(ir),'ko')
%%