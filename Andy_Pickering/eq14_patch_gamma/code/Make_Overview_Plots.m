%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_Overview_Plots.m
%
% Making plots for chi-pod overview notes
%
% Depends:
% - Get_and_bin_profiles.m
% - chi_vs_eps_normalized_plot.m
% - ComputeGamma.m
% - get_cham_cnums_eq14.m
%
%
%---------------
% 4/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;
dz = 2 ;
cnums_to_get = get_cham_cnums_eq14;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

screen_chi=1

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,screen_chi);


%% Pcolor of chipod & cham chi, and N2,Tz

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hline(80,'k--')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hline(80,'k--')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])

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
hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chameleon')
ylabel('P [db]')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chi-pod')
ylabel('P [db]')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hline(80,'k--')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hline(80,'k--')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])

%figname = [project_short '_Pcolor_BothEps_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
figname = [project_short '_Pcolor_BothEps_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir, figname), '-dpng')

%% Plot chipod method vs chameleon

icham = find(cham.P>80);
ichi = find(chipod.P>80);

figure(3);clf
agutwocolumn(1)
wysiwyg

subplot(211)
histogram2( log10(cham.chi(icham,:)), log10(chipod.chi(ichi,:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-12 -4])
ylim([-12 -4])
xlabel('\chi chameleon')
ylabel('\chi chipod')

subplot(212)
histogram2( log10(cham.eps(icham,:)), log10(chipod.eps(ichi,:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')

xlim([-8.5 -4])
ylim([-8.5 -4])
xlabel('\epsilon chameleon')
ylabel('\epsilon chipod')

%figname = [project_short '_chamVschipod_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
figname = [project_short '_chamVschipod_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir, figname), '-dpng')


%% Pcolor of *ratio* of chipod & cham eps, plus N2,Tz

figure(4);clf
agutwocolumn(1)
wysiwyg

rr=2 ;
cc=1 ;

ax1 = subplot(rr,cc,1) ;
ezpc(chipod.cnum,chipod.P,(chipod.chi ./ cham.chi))
caxis([0.1 1.5])
colorbar
ylabel('P [db]')
hline(80,'k--')

ax2 = subplot(rr,cc,2);
ezpc(cham.cnum,cham.P,(chipod.eps ./ cham.eps))
caxis([0 1.5])
colorbar
%title('\epsilon chi-pod')
ylabel('P [db]')
hline(80,'k--')

linkaxes([ax1 ax2 ])

%figname = [project_short '_Pcolor_BothEps_N2_Tz_zsmooth_' num2str(Params.z_smooth) '_' num2str(dz) 'mbin']
%print(fullfile(fig_dir, figname), '-dpng')


%% chi vs eps, normalized by N2/Tz etc

figure(5);clf
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


%% See if gamma computed from multi-profile averageds of N2,Tz,chi,eps is 0.2?
% compare to gamma computed from individual 1m data points in every profile

clear ; close all

dz=20; % bin size to average over

eq14_patches_paths

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

cnum_range = [2400 3000];

clear cnums
cnums = [cnum_range(1) : cnum_range(2) ];

bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums = setdiff(cnums,bad_prof);

iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));
%
eps = cham.EPSILON(:,iCham); eps = eps(:) ;
chi = cham.CHI(:,iCham)  ; chi = chi(:) ;
N2  = cham.N2(:,iCham)   ; N2 = N2(:) ;
Tz  = cham.DTDZ(:,iCham) ; Tz = Tz(:) ;
P   = cham.P(:,iCham)    ; P = P(:) ;

% compute gamma from these values
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
gam = ComputeGamma(N2,Tz,chi,eps);

% Now average profiles together in 10m bins and compute gamma from that

[eps_bin zout Nobs] = binprofile(eps ,P, 0, dz, 200,1);
[chi_bin zout Nobs] = binprofile(chi ,P, 0, dz, 200,1);
[N2_bin  zout Nobs] = binprofile(N2  ,P, 0, dz, 200,1);
[Tz_bin  zout Nobs] = binprofile(Tz  ,P, 0, dz, 200,1);

gam_avg = ComputeGamma(N2_bin,Tz_bin,chi_bin,eps_bin);

figure(6);clf
agutwocolumn(0.6)
wysiwyg

ax1 = subplot(121) ;
boxplot(log10(gam))
hline(log10(0.2),'k--')
grid on
ylabel('log_{10}[\gamma]','fontsize',16)
title(['1mavg, profiles ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])

ax2 = subplot(122) ;
boxplot(log10(gam_avg))
grid on
hline(log10(0.2),'k--')
title(['profile-averaged, ' num2str(dz) ' m binned'])

linkaxes([ax1 ax2])

figname=[project_short '_gamma_point_avg_box_' num2str(dz) 'mbinned']
print(fullfile(fig_dir,figname),'-dpng')


%%

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7   ;
Params.z_smooth = 10 ;

screen_chi=1

%dz = 10 % bin size
zmin=0  ;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths
Pmin = 80;
cnums_to_get = get_cham_cnums_eq14 ;
%cnums_to_get = 2000:3000;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

figure(7);clf
agutwocolumn(1)
wysiwyg

iax=1
rr=3
cc=2
for dz=[1 10 50]
    
    clear chipod cham
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,screen_chi)
    
    % Nan values in mixed layer
    clear ib
    ib = find(cham.P<Pmin);
    cham.eps(ib) = nan;
    
    clear ib
    ib = find(chipod.P<Pmin);
    chipod.eps(ib) = nan;
    
    subplot(rr,cc,iax)
    h = chi_vs_eps_normalized_plot(cham.eps, cham.chi, cham.N2, cham.Tz)
    title([project_short ' Chameleon ' num2str(dz) 'm binned, >' num2str(Pmin) 'db'])
    
    iax=iax+1;
    
    subplot(rr,cc,iax)
    hh=histogram2(  real(log10(cham.eps)),log10(chipod.eps),80,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec-1,'r--')
    plot(xvec,xvec+1,'r--')
    
    if screen_chi==1
        ylim([-8.5 -4])
        xlim([-8.5 -4])
    else
        ylim([-11 -4])
        xlim([-11 -4])
    end
    ylabel('log_{10} [\epsilon_{\chi}]','fontsize',16)
    if iax>6
        xlabel('log_{10} [\epsilon ]','fontsize',16)
    end
    
    iax=iax+1;
    
end

%
figname=['eq14_NormScat_chiVscham_diff_dz_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir,figname),'-dpng')

%% plot chi vs chi and eps vs eps for different depth bin averaging

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7   ;
Params.z_smooth = 10 ;

screen_chi=1

%dz = 10 % bin size
zmin=0  ;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths
Pmin = 80;
cnums_to_get = get_cham_cnums_eq14 ;
%cnums_to_get = 2000:3000;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

figure(7);clf
agutwocolumn(1)
wysiwyg

iax=1
rr=3
cc=2
for dz=[1 10 50]
    
    clear chipod cham
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,screen_chi)
    
    % Nan values in mixed layer
    clear ib
    ib = find(cham.P<Pmin);
    cham.eps(ib) = nan;
    
    clear ib
    ib = find(chipod.P<Pmin);
    chipod.eps(ib) = nan;
    
    subplot(rr,cc,iax)
    %    h = chi_vs_eps_normalized_plot(cham.eps, cham.chi, cham.N2, cham.Tz)
    hh=histogram2(  real(log10(cham.chi)),log10(chipod.chi),80,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec-1,'r--')
    plot(xvec,xvec+1,'r--')
    ylim([-10 -4])
    xlim([-10 -4])
    ylabel('log_{10} [\chi_{\chi}]','fontsize',16)
    
    title([project_short ' Chameleon ' num2str(dz) 'm binned, >' num2str(Pmin) 'db'])
    if iax==rr*2 -1
        xlabel('log_{10} [\chi ]','fontsize',16)
    end
    
    iax=iax+1;
    
    subplot(rr,cc,iax)
    hh=histogram2(  real(log10(cham.eps)),log10(chipod.eps),80,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec-1,'r--')
    plot(xvec,xvec+1,'r--')
    
    if screen_chi==1
        ylim([-8.5 -4])
        xlim([-8.5 -4])
    else
        ylim([-11 -4])
        xlim([-11 -4])
    end
    ylabel('log_{10} [\epsilon_{\chi}]','fontsize',16)
    if iax==rr*2
        xlabel('log_{10} [\epsilon ]','fontsize',16)
    end
    
    iax=iax+1;
    
end

%
figname=['eq14_chiVscham_chiANDeps_diff_dz_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir,figname),'-dpng')


%% histogram of epsilon ratio for different size depth averaging

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7   ;
Params.z_smooth = 10 ;

screen_chi=1

%dz = 10 % bin size
zmin=0  ;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths
Pmin = 80;
cnums_to_get = get_cham_cnums_eq14 ;
%cnums_to_get = 2000:3000;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

figure(7);clf
agutwocolumn(1)
wysiwyg

iax=1
cols=['b','r','y','m']
h=[]

for dz = [1 10 50]
    
    clear chipod cham
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,screen_chi)
    
    % Nan values in mixed layer
    clear ib
    ib = find(cham.P<Pmin);
    cham.eps(ib) = nan;
    
    clear ib
    ib = find(chipod.P<Pmin);
    chipod.eps(ib) = nan;

    subplot(2,1,1)
    hh=histogram( log10( chipod.chi(:) ./ cham.chi(:)),[-3:0.1:3], 'Normalization','pdf','Edgecolor','none','FaceAlpha',0.5)
    grid on
    xlim([-2 2])
    ylim([0 1.3])
    hold on
    freqline(nanmean(log10( chipod.chi(:) ./ cham.chi(:))),cols(iax))
    hold on
    h=[h hh];
%    iax=iax+1;

    subplot(2,1,2)
    hh=histogram( log10( chipod.eps(:) ./ cham.eps(:)),[-3:0.1:3], 'Normalization','pdf','Edgecolor','none','FaceAlpha',0.5)
    grid on
    xlim([-3 2])
    ylim([0 0.9])
    hold on
    freqline(nanmean(log10( chipod.eps(:) ./ cham.eps(:))),cols(iax))
    hold on
    h=[h hh];
    iax=iax+1;
    
end


legend(h,'1m','10m','50m')
xlabel(['\epsilon_{\chi}/\epsilon'],'fontsize',16)
%
subplot(211)
xlabel(['\chi_{\chi}/\chi'],'fontsize',16)

%
figname=['eq14_chiVscham_hist_diff_dz_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir,figname),'-dpng')


%% make similar plot, but for averaging diffferent numbers of profiles

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth =10 ;

screen_chi=1

dz = 10 % bin size
Pmin = 80

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(8);clf
agutwocolumn(1)
wysiwyg

iax=1;
for dp = [2 10 50]
    
    eps_cham_all = [];
    chi_cham_all = [];
    N2_cham_all = [];
    Tz_cham_all = [];
    
    eps_chi_all = [];
    chi_chi_all = [];
    N2_chi_all = [];
    Tz_chi_all = [];
    
    
    for ix = 1:round(3000/dp)%
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get = setdiff(cnums_to_get,bad_prof);
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi);
        
        eps_cham_all = [eps_cham_all(:) ; eps_cham_avg(:) ] ;
        chi_cham_all = [chi_cham_all(:) ; chi_cham_avg(:) ] ;
        N2_cham_all  = [N2_cham_all(:)  ; N2_cham_avg(:) ] ;
        Tz_cham_all  = [Tz_cham_all(:)  ; Tz_cham_avg(:) ] ;
        
        eps_chi_all = [eps_chi_all(:) ; eps_chi_avg(:) ] ;
        chi_chi_all = [chi_chi_all(:) ; chi_chi_avg(:) ] ;
        N2_chi_all  = [N2_chi_all(:)  ; N2_chi_avg(:) ] ;
        Tz_chi_all  = [Tz_chi_all(:)  ; Tz_chi_avg(:) ] ;
        
        
    end % idx
    
    subplot(3,2,iax)
    h = chi_vs_eps_normalized_plot(eps_cham_all, chi_cham_all, N2_cham_all, Tz_cham_all)
    title([num2str(dp) ' profile averages'])
    
    iax = iax+1;
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(eps_cham_all)),log10(eps_chi_all),20,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    
    if screen_chi==1
        ylim([-8.5 -5]); xlim([-8.5 -5])
    else
        ylim([-11 -4]); xlim([-11 -4])
    end
    
    ylabel('log_{10} [\epsilon_{\chi}]','fontsize',16)
    
    if iax>4
        xlabel('log_{10} [\epsilon ]','fontsize',16)
    end
    
    title([num2str(dp) ' profile averages'])
    
    iax = iax+1;
    
end % dp

%
figname=['eq14_NormScat_chiVscham_diff_prof_avg_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir,figname),'-dpng')

%% Plot chi vs chi, eps vs eps, for different # profiles averaged

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth =10 ;

screen_chi=1

dz = 10 % bin size
Pmin = 80

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(8);clf
agutwocolumn(1)
wysiwyg

iax=1;
for dp = [2 10 50]
    
    eps_cham_all = [];
    chi_cham_all = [];
    N2_cham_all = [];
    Tz_cham_all = [];
    
    eps_chi_all = [];
    chi_chi_all = [];
    N2_chi_all = [];
    Tz_chi_all = [];
    
    
    for ix = 1:round(3000/dp)%
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get = setdiff(cnums_to_get,bad_prof);
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi);
        
        eps_cham_all = [eps_cham_all(:) ; eps_cham_avg(:) ] ;
        chi_cham_all = [chi_cham_all(:) ; chi_cham_avg(:) ] ;
        N2_cham_all  = [N2_cham_all(:)  ; N2_cham_avg(:)  ] ;
        Tz_cham_all  = [Tz_cham_all(:)  ; Tz_cham_avg(:)  ] ;
        
        eps_chi_all = [eps_chi_all(:) ; eps_chi_avg(:) ] ;
        chi_chi_all = [chi_chi_all(:) ; chi_chi_avg(:) ] ;
        N2_chi_all  = [N2_chi_all(:)  ; N2_chi_avg(:)  ] ;
        Tz_chi_all  = [Tz_chi_all(:)  ; Tz_chi_avg(:)  ] ;
        
        
    end % idx
    
    subplot(3,2,iax)
    %    h = chi_vs_eps_normalized_plot(eps_cham_all, chi_cham_all, N2_cham_all, Tz_cham_all)
    hh=histogram2(  real(log10(chi_cham_all)),log10(chi_chi_all),40,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    ylim([-10 -4]); xlim([-10 -4])
    
    title([num2str(dp) ' profile averages'])
    ylabel('log_{10} [\chi_{\chi}]','fontsize',16)
    
    if iax==5
        xlabel('log_{10} [\chi ]','fontsize',16)
    end
    
    
    iax = iax+1;
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(eps_cham_all)),log10(eps_chi_all),25,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    
    if screen_chi==1
        ylim([-8.5 -5]); xlim([-8.5 -5])
    else
        ylim([-11 -4]); xlim([-11 -4])
    end
    
    ylabel('log_{10} [\epsilon_{\chi}]','fontsize',16)
    
    
    if iax==6
        xlabel('log_{10} [\epsilon ]','fontsize',16)
    end
    
    title([num2str(dp) ' profile averages'])
    
    iax = iax+1;
    
end % dp

%
figname=['eq14_chiVscham_chiANDeps_diff_prof_avg_screen_chi_' num2str(screen_chi)]
print(fullfile(fig_dir,figname),'-dpng')


%% May 1 2017 - plot histograms of eps_chi/eps for differnt # prof avg.

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth =10 ;

screen_chi=1

dz = 10 % bin size
Pmin = 80

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(8);clf
agutwocolumn(1)
wysiwyg

h=[]
iax=1
cols=['b','r','m']
for dp = [1 10 50 ]
    
    eps_cham_all = [];
    chi_cham_all = [];
    N2_cham_all  = [];
    Tz_cham_all  = [];
    
    eps_chi_all = [];
    chi_chi_all = [];
    N2_chi_all  = [];
    Tz_chi_all  = [];
    
    for ix = 1:round(3000/dp)%
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get = setdiff(cnums_to_get,bad_prof);
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi);
        
        eps_cham_all = [eps_cham_all(:) ; eps_cham_avg(:) ] ;
        chi_cham_all = [chi_cham_all(:) ; chi_cham_avg(:) ] ;
        N2_cham_all  = [N2_cham_all(:)  ; N2_cham_avg(:) ] ;
        Tz_cham_all  = [Tz_cham_all(:)  ; Tz_cham_avg(:) ] ;
        
        eps_chi_all = [eps_chi_all(:) ; eps_chi_avg(:) ] ;
        chi_chi_all = [chi_chi_all(:) ; chi_chi_avg(:) ] ;
        N2_chi_all  = [N2_chi_all(:)  ; N2_chi_avg(:) ] ;
        Tz_chi_all  = [Tz_chi_all(:)  ; Tz_chi_avg(:) ] ;
        
    end % idx

    subplot(211)
    histogram(  log10(chi_chi_all./chi_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    xlim([-2 2])
    ylim([0 1.5])
    grid on
    hold on
    freqline(nanmean(log10(chi_chi_all./chi_cham_all)),cols(iax))
    text(1.5,0.8-(iax*0.1),num2str(roundx(nanmean(log10(chi_chi_all./chi_cham_all)),2)),'color',cols(iax),'fontsize',14)
    hold on

    subplot(212)
    hh=histogram(  log10(eps_chi_all./eps_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax))
    xlim([-2 2])
    ylim([0 1.2])
    grid on
    hold on
    freqline(nanmean(log10(eps_chi_all./eps_cham_all)),cols(iax))
    text(1.5,0.8-(iax*0.1),num2str(roundx(nanmean(log10(eps_chi_all./eps_cham_all)),2)),'color',cols(iax),'fontsize',14)
    hold on
    
    h=[h hh];
    iax=iax+1
    
end % dp
%%
subplot(212)
legend(h,'1 profile','10 profiles','50 profiles')
xlabel('\epsilon_{\chi}/\epsilon','fontsize',16)
ylabel('pdf','fontsize',16)

subplot(211)
xlabel('\chi_{\chi}/\chi','fontsize',16)
ylabel('pdf','fontsize',16)
%

print( fullfile(fig_dir,'eq14_eps_ratio_hist_diff_prof_avg'), '-dpng')

%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

screen_chi = 1

dz = 10 % bin size
zmin=0;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths

cnums_to_get = get_cham_cnums_eq14 ;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,screen_chi)

% Nan values in mixed layer

Pmin = 80;

clear ib
ib = find(cham.P<Pmin);
cham.eps(ib) = nan;

clear ib
ib = find(chipod.P<Pmin);
chipod.eps(ib) = nan;

%
% Histograms of ratio of chi-pod epsilon to chameleon epsilon (10bins)

figure(9);clf
h1 = histogram(log10( chipod.eps(:) ./ cham.eps(:) ),'EdgeColor','none','Normalization','pdf');
hold on
xlim([-3 3])
grid on
xlabel('log_{10}[\epsilon_{\chi} /\epsilon ]')
ylabel('pdf')
title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])
freqline(nanmean(log10( chipod.eps(:) ./ cham.eps(:) )))

%figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_Pmin' num2str(Pmin)]figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_Pmin' num2str(Pmin)]
figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_Pmin' num2str(Pmin) '_screen_chi_' num2str(screen_chi)]
print( fullfile(fig_dir,figname),'-dpng')

%

% Plot eps vs chi, normalized so slope is equal to 1/2*gamma

figure(10);clf
agutwocolumn(0.8)
wysiwyg

hh=histogram2(  real(log10(cham.eps./cham.N2)),log10(cham.chi./(cham.Tz.^2)),80,'DisplayStyle','tile')
grid on
hold on
xvec=linspace(1e-7,1e-1,100);
h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
ylim([-7.5 -1])
xlim([-5.5 -1])
ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
title([project_short ' Chameleon ' num2str(dz) 'm binned, >' num2str(Pmin) 'db'])

fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_Pmin_' num2str(Pmin)]
print( fullfile(fig_dir,fname),'-dpng')

%%