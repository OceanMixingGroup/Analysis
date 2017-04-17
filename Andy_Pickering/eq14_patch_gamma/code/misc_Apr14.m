%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
% Making plots for chi-pod paper outline
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;
dz = 2 ;
cnums_to_get = get_cham_cnums_eq14;

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200);

%% Pcolor of chipod & cham eps, and N2,Tz

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1) ;
ezpc(cham.cnum,cham.P,log10(cham.eps))
caxis([-11 -4])
colorbar
title('\epsilon chameleon')
ylabel('P [db]')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
caxis([-11 -4])
colorbar
title('\epsilon chi-pod')
ylabel('P [db]')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('dT/dz')

linkaxes([ax1 ax2 ax3 ax4])

figname = [project_short '_Pcolor_BothEps_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
print(fullfile(fig_dir, figname), '-dpng')

%% Pcolor of chipod & cham chi, and N2,Tz

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
caxis([-11 -4])
colorbar
title('\chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
caxis([-11 -4])
colorbar
title('\chi chameleon')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('dT/dz')

linkaxes([ax1 ax2 ax3 ax4])

figname = [project_short '_Pcolor_BothChi_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
print(fullfile(fig_dir, figname), '-dpng')


%% 2D histogrmas of chipod vs chameleon

figure(2);clf
agutwocolumn(1)
wysiwyg

subplot(211)
histogram2( log10(cham.chi(:)), log10(chipod.chi(:)),'DisplayStyle','tile')
grid on
xlim([-12 -3])
ylim([-12 -3])
hold on
xvec = linspace(-12,-3,100);
plot(xvec,xvec,'k--')
xlabel('\chi chameleon','fontsize',16)
ylabel('\chi chi-pod','fontsize',16)

subplot(212)
histogram2( log10(cham.eps(:)), log10(chipod.eps(:)),80,'DisplayStyle','tile')
grid on
xlim([-11 -4])
ylim([-11 -4])
hold on
xvec = linspace(-12,-3,100);
plot(xvec,xvec,'k--')
xlabel('\epsilon chameleon','fontsize',16)
ylabel('\epsilon chi-pod','fontsize',16)

%% chi vs eps, normalized by N2/Tz etc

figure(2);clf
agutwocolumn(0.8)
wysiwyg

Pmin=80;
iz = find(cham.P>Pmin);

hh=histogram2(  real(log10(cham.eps(iz,:)./cham.N2(iz,:))),log10(cham.chi(iz,:)./(cham.Tz(iz,:).^2)),80,'DisplayStyle','tile')
grid on
hold on
xvec=linspace(1e-7,1e-1,100);
h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
ylim([-7 -1])
xlim([-7 -1])
ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
title([project_short ' Chameleon 10m binned, >' num2str(Pmin) 'db'])


%% Now do same plot, but average some profiles together?

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;

dz = 10 % bin size
Pmin = 0

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

cnums_to_get = 2500:3000;

[eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
    Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short)

%[eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
%    Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short)
%

figure(1);clf
hcham = plot(log10(eps_cham_avg),P_cham) ;
hold on
hchi = plot(log10(eps_chi_avg),P_chi);
axis ij
xlim([-10 -4])
legend([hcham hchi],'cham','chi')
grid on
%%
