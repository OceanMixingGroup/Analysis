%
%
%
%
%
%
%%
clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32  ;
Params.z_smooth=10;
dz = 2 ;

cnums_to_get = 200:1000;%get_cham_cnums_eq14;
%bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
%cnums_to_get = setdiff(cnums_to_get,bad_prof);

eq08_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

screen_chi=1
screen_ml=0
Pmin=0;

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

%load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat')

%
%% Pcolor of chipod & cham chi, and N2,Tz

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])
%
figname = [project_short '_Pcolor_BothChi_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir, figname), '-dpng')


%% Pcolor of chipod & cham eps, and N2,Tz

figure(2);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1) ;
ezpc(cham.cnum,cham.P,log10(cham.eps))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chameleon')
ylabel('P [db]')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chi-pod')
ylabel('P [db]')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])
%
figname = [project_short '_Pcolor_BothEps_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir, figname), '-dpng')



%% Plot 2D histogram of chipod method vs chameleon

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32  ;
Params.z_smooth=10;

dz = 2 ;

cnums_to_get = 200:1000%get_cham_cnums_eq14;
%bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
%cnums_to_get = setdiff(cnums_to_get,bad_prof);

eq08_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

screen_chi=1
Pmin=80
screen_ml=0

% reload data, screening convective regions
[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat')

figure(3);clf
agutwocolumn(1)
wysiwyg

subplot(211)
%histogram2( log10(cham.chi(icham,:)), log10(chipod.chi(ichi,:)), 'DisplayStyle','tile')
histogram2( log10(cham.chi(:)), log10(chipod.chi(:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-12 -4])
ylim([-12 -4])
xlabel('\chi','fontsize',16)
ylabel('\chi_{\chi}','fontsize',16)

subplot(212)
%histogram2( log10(cham.eps(icham,:)), log10(chipod.eps(ichi,:)),50, 'DisplayStyle','tile')
histogram2( log10(cham.eps(:)), log10(chipod.eps(:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-8.5 -4])
ylim([-8.5 -4])
xlabel('\epsilon ','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)

%
%figname = [project_short '_chamVschipod_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
figname = [project_short '_chamVschipod_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir, figname), '-dpng')


%%