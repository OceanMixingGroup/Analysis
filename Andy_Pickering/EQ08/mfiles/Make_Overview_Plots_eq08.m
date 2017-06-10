%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_Overview_Plots_eq08.m
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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma    = 0.2;
Params.fmax     = 10  ;
Params.z_smooth = 10 ;
Params.resp_corr= 0  ;
Params.fc       = 99 ;

eq08_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/
%
dz = 2 ;
cnums_to_get = 200:2700;

screen_chi= 1 ;
screen_ml = 0 ;
Pmin      = 0 ;

%project = 'eq08'
%eval([project '_patches_paths'])

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
if exist(fullfile(analysis_dir,project_short,'Data',save_name),'file')==2
load(fullfile(analysis_dir,project_short,'Data',save_name))
else
[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq08_patch_gamma/data/eq08_mldepths.mat')


%% Pcolor of chipod & cham chi, and N2,Tz

h = pcolor_chi(chipod,cham)  ;
%
for i=1:4
    subplot(4,1,i)
    hold on
    plot(zml_cnum,zml,'k')
end
shg

figname = [project_short '_Pcolor_BothChi_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]

print(fullfile(fig_dir, figname), '-dpng')

%% Pcolor of chipod & cham *eps*, and N2,Tz

h = pcolor_eps(chipod,cham) ;

%
for i=1:4
    subplot(4,1,i)
    hold on
    plot(zml_cnum,zml,'k')
end
shg

figname = [project_short '_Pcolor_BothEps_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]

print(fullfile(fig_dir, figname), '-dpng')

%

%% Combine above 2 plots - plot chi,eps,n2,dtdz on 6X1 panel plot

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=6 ; cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3) ;
ezpc(cham.cnum,cham.P,log10(cham.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chameleon')
ylabel('P [db]')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chi-pod')
ylabel('P [db]')

ax5 = subplot(rr,cc,5);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
plot(zml_cnum,zml,'k')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax6 = subplot(rr,cc,6);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
plot(zml_cnum,zml,'k')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6])
figname = [project_short '_Pcolor_Both_epsANDChi_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)] ;
print(fullfile(fig_dir, figname), '-dpng')


%% 2D hist vs depth?

%clear cham chipod ;

% Params.gamma    = 0.2;
% Params.fmax     = 15 ;
% Params.z_smooth = 10 ;
% Params.resp_corr= 0  ;
% Params.fc       = 99 ;
%
% dz = 2 ;

%cnums_to_get = 200:2800 ;

% screen_chi = 1 ;
% Pmin       = 0 ;
% screen_ml  = 0 ;

eq08_patches_paths

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

%[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

P = repmat(chipod.P,1,length(chipod.cnum));

figure(33);clf
agutwocolumn(1)
wysiwyg

subplot(211)
histogram2( real(log10( chipod.chi ./ cham.chi )), P, 'DisplayStyle','tile','EdgeColor','none')
axis ij
ylabel('P[db]','fontsize',16)
xlabel('log_{10}[\chi_{\chi}/\chi]','fontsize',16)
xlim([-3 3])
ylim([0 200])
colorbar
caxis([0 1000])
freqline(0,'k-')
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])

subplot(212)
histogram2( real(log10( chipod.eps ./ cham.eps )), P, 'DisplayStyle','tile','EdgeColor','none')
axis ij
ylabel('P[db]','fontsize',16)
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]','fontsize',16)
xlim([-3 3])
ylim([0 200])
caxis([0 1000])
colorbar
freqline(0,'k-')
%
figname = [project_short '_chi_eps_Vs_P_2Dhist_screen_chi_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')


%% Plot 2D histograms of chipod method vs chameleon, for chi and eps

clear cham chipod

% Params.gamma     = 0.2;
% Params.fmax      = 15 ;
% Params.z_smooth  = 10 ;
% Params.resp_corr = 0  ;
% Params.fc        = 99 ;
%
% dz = 2 ;

% Reload data and get rid of mixed layer regions
cnums_to_get = 200:2700;
screen_chi = 1 ;
Pmin       = 20;
screen_ml  = 1 ;
dz = 2 ;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
if exist(fullfile(analysis_dir,project_short,'Data',save_name))==2
load(fullfile(analysis_dir,project_short,'Data',save_name))
else
[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end

%h = scatter_chi_eps_chipod_cham(chipod,cham) ;

h = figure;clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(211) ;
histogram2( log10(cham.chi(:)), log10(chipod.chi(:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-11 -4])
ylim([-11 -4])
xlabel('log_{10}[\chi]','fontsize',16)
ylabel('log_{10}[\chi_{\chi}]','fontsize',16)
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])

ax2 = subplot(212);
histogram2( log10(cham.eps(:)), log10(chipod.eps(:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-8.5 -4.5])
ylim([-8.5 -4.5])
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)

figname = [project_short '_chamVschipod_screen_chi_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir, figname), '-dpng')

%%
%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

%clear cham chipod;% close all

% Params.gamma    = 0.2 ;
% Params.fmax     = 15  ;
% Params.z_smooth = 10  ;
% Params.resp_corr= 0   ;     % correct TP spectra for freq response of thermistor?
% Params.fc       = 99  ;    % cutoff frequency for response correction

% Pmin       = 20
% screen_ml  = 1
% screen_chi = 1

%dz = 10 % bin size

%zmin=0;
%zmax=200;

%addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
%addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

%eq08_patches_paths

%cnums_to_get = 200:2800;
%[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,Pmin,screen_chi,screen_ml) ;

% Histograms of ratio of chi-pod epsilon to chameleon epsilon (10bins)

figure(13);clf
agutwocolumn(0.8)
wysiwyg

subplot(211)

h1 = histogram(log10( chipod.chi(:) ./ cham.chi(:) ),'EdgeColor','none','Normalization','pdf') ;
hold on
xlim([-2 2])
grid on
xlabel('log_{10}[\chi_{\chi} /\chi ]')
ylabel('pdf')
%title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])
%title([project_short ' ' num2str(dz) ' m binned '])
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])
freqline( nanmean(log10( chipod.chi(:) ./ cham.chi(:) )))
text(1,0.8,['\mu= ' num2str( roundx(nanmean(log10( chipod.chi(:) ./ cham.chi(:) )),2))],'fontsize',15)

subplot(212)
h2 = histogram(log10( chipod.eps(:) ./ cham.eps(:) ),'EdgeColor','none','Normalization','pdf') ;
hold on
xlim([-3 3])
grid on
xlabel('log_{10}[\epsilon_{\chi} /\epsilon ]')
ylabel('pdf')
%title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])
%title([project_short ' ' num2str(dz) ' m binned '])
freqline( nanmean(log10( chipod.eps(:) ./ cham.eps(:) )))
text(1,0.4,['\mu= ' num2str( roundx(nanmean(log10( chipod.eps(:) ./ cham.eps(:) )),2))],'fontsize',15)


clear figname
figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_screen_chi_' num2str(screen_chi) '_screenml_' num2str(screen_ml)  '_' MakeChiPathStr(Params)]
print( fullfile(fig_dir,figname),'-dpng')
clear figname


%% Plot eps vs chi, normalized so slope is equal to 1/2*gamma

figure(14);clf
agutwocolumn(0.8)
wysiwyg

hh = histogram2(  real(log10(cham.eps./cham.N2)),log10(cham.chi./(cham.Tz.^2)),80,'DisplayStyle','tile') ;
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
%title([project_short ' Chameleon ' num2str(dz) 'm binned, >' num2str(Pmin) 'db'])
%title([project_short ' Chameleon ' num2str(dz) 'm binned'])
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])

fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_' MakeChiPathStr(Params)]
print( fullfile(fig_dir,fname),'-dpng')



%% plot chi vs chi and eps vs eps for different depth bin averaging

clear cham chipod; %close all

addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

% Params.gamma    = 0.2 ;
% Params.fmax     = 15  ;
% Params.z_smooth = 10   ;
% Params.resp_corr= 0   ;  % correct TP spectra for freq response of thermistor?
% Params.fc       = 99  ;  % cutoff frequency for response correction

screen_chi = 1  ;
screen_ml  = 1  ;
Pmin       = 20 ;

zmin=0  ;
zmax=200;

cnums_to_get = 200:2700;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq08_patches_paths

figure(8);clf
agutwocolumn(1)
wysiwyg

iax = 1 ;
rr=3 ; cc=2 ;

for dz=[2 10 50]
    
    save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
    if exist(fullfile(analysis_dir,project_short,'Data',save_name),'file')==2
        load(fullfile(analysis_dir,project_short,'Data',save_name))
    else
        [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
    end
    
    subplot(rr,cc,iax)
    hh = histogram2(  real(log10(cham.chi)),log10(chipod.chi),'XBinEdges',[-10:0.15:-4],'YBinEdges',[-10:0.15:-4],'DisplayStyle','tile','EdgeColor','none')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec-1,'r--')
    plot(xvec,xvec+1,'r--')
    ylim([-10 -4])
    xlim([-10 -4])
    ylabel('log_{10} [\chi_{\chi}]','fontsize',16)
    
    title([project_short ' ' num2str(dz) 'm binned '])
    if iax==rr*2 -1
        xlabel('log_{10} [\chi ]','fontsize',16)
    end
    
    iax=iax+1;
    
    subplot(rr,cc,iax)
    hh = histogram2(  real(log10(cham.eps)),log10(chipod.eps),'XBinEdges',[-8.5:0.15:-4],'YBinEdges',[-8.5:0.15:-4],'DisplayStyle','tile','EdgeColor','none') ;
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
clear figname
figname=[project_short '_chiVscham_chiANDeps_diff_dz_screen_chi_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')
clear figname

%% histogram of epsilon ratio for different size depth averaging

clear cham chipod; %close all

% Params.gamma    = 0.2 ;
% Params.fmax     = 15  ;
% Params.z_smooth = 10  ;
% Params.resp_corr= 0   ;  % correct TP spectra for freq response of thermistor?
% Params.fc       = 99  ;  % cutoff frequency for response correction

screen_chi = 1 ;
Pmin       = 20 ;
screen_ml  = 1 ;

zmin=0  ;
zmax=200;

cnums_to_get = 200:2700;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq08_patches_paths


figure(9);clf
agutwocolumn(1)
wysiwyg

iax  = 1 ;
cols = ['b','r','g'];
h = [] ;

for dz = [1 10 50]
    
    clear chipod cham
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,Pmin,screen_chi,screen_ml) ;
    
    subplot(2,1,1)
    hh=histogram( log10( chipod.chi(:) ./ cham.chi(:)),[-3:0.1:3], 'Normalization','pdf','DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    grid on
    xlim([-2 2])
    %ylim([0 1.4])
    hold on
    freqline(nanmean(log10( chipod.chi(:) ./ cham.chi(:))),cols(iax))
    hold on
    chirat=nanmean(log10( chipod.chi(:) ./ cham.chi(:)));
    text(1,0.8-(iax*0.1),['\mu= ' num2str(roundx(chirat,2))],'fontsize',15,'color',cols(iax))

    h=[h hh];
    
    subplot(2,1,2)
    histogram( log10( chipod.eps(:) ./ cham.eps(:)),[-3:0.1:3], 'Normalization','pdf','DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    grid on
    xlim([-3 3])
    ylim([0 1])
    hold on
    freqline(nanmean(log10( chipod.eps(:) ./ cham.eps(:))),cols(iax))
    hold on
    %h=[h hh];
    epsrat=nanmean(log10( chipod.eps(:) ./ cham.eps(:)));
    text(1,0.8-(iax*0.1),['\mu= ' num2str(roundx(epsrat,2))],'fontsize',15,'color',cols(iax))
    
    iax=iax+1;
    
end


legend(h,'1m','10m','50m')
xlabel(['\epsilon_{\chi}/\epsilon'],'fontsize',16)
ylabel('pdf','fontsize',16)
%
subplot(211)
xlabel(['\chi_{\chi}/\chi'],'fontsize',16)
ylabel('pdf','fontsize',16)
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ])

%
clear figname
figname=[project_short '_chiVscham_hist_diff_dz_screen_chi_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')
clear figname


%% Plot chi vs chi, eps vs eps, for different # profiles averaged

clear cham chipod ; %close all

% Params.gamma    = 0.2 ;
% Params.fmax     = 15  ;
% Params.z_smooth = 10  ;
% Params.resp_corr= 0   ; % correct TP spectra for freq response of thermistor?
% Params.fc       = 99  ; % cutoff frequency for response correction


screen_chi= 1
screen_ml = 1
Pmin      = 20

dz        = 10 % bin size

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths

figure(11);clf
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
    N2_chi_all  = [];
    Tz_chi_all  = [];
    
    
    for ix = 1:round(2800/dp)%
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [chipod, cham] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
        
        eps_cham_all = [eps_cham_all(:) ; cham.eps(:) ] ;
        chi_cham_all = [chi_cham_all(:) ; cham.chi(:) ] ;
        N2_cham_all  = [N2_cham_all(:)  ; cham.N2(:) ] ;
        Tz_cham_all  = [Tz_cham_all(:)  ; cham.Tz(:) ] ;
        
        eps_chi_all = [eps_chi_all(:) ; chipod.eps(:) ] ;
        chi_chi_all = [chi_chi_all(:) ; chipod.chi(:) ] ;
        N2_chi_all  = [N2_chi_all(:)  ; chipod.N2(:) ] ;
        Tz_chi_all  = [Tz_chi_all(:)  ; chipod.Tz(:) ] ;
        
    end % idx
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(chi_cham_all)),log10(chi_chi_all),40,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec+1,'r--')
    plot(xvec,xvec-1,'r--')
    
    ylim([-10 -4]); xlim([-10 -4])
    
    title([num2str(dp) ' profile averages'])
    ylabel('log_{10} [\chi_{\chi}]','fontsize',16)
    
    if iax==5
        xlabel('log_{10} [\chi ]','fontsize',16)
    end
    
    
    iax = iax+1;
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(eps_cham_all)),log10(eps_chi_all),20,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    plot(xvec,xvec+1,'r--')
    plot(xvec,xvec-1,'r--')
    
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
clear figname
figname=[project_short '_chiVscham_chiANDeps_diff_prof_avg_screen_chi_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')


%% May 1 2017 - plot histograms of eps_chi/eps for differnt # prof avg.

clear cham chipod;% close all
%
% Params.gamma    = 0.2 ;
% Params.fmax     = 15  ;
% Params.z_smooth = 10  ;
% Params.resp_corr= 0   ;  % correct TP spectra for freq response of thermistor?
% Params.fc       = 99  ;  % cutoff frequency for response correction

screen_chi = 1  ;
Pmin       = 20 ;
screen_ml  = 1  ;

dz = 10 % bin size

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths

figure(12);clf
agutwocolumn(1)
wysiwyg

h=[];
iax=1;
cols=['b','r','g'];
for dp = [1 10 50 ]
    
    eps_cham_all = [];
    chi_cham_all = [];
    N2_cham_all  = [];
    Tz_cham_all  = [];
    
    eps_chi_all = [];
    chi_chi_all = [];
    N2_chi_all  = [];
    Tz_chi_all  = [];
    
    for ix = 1 : round(2800/dp)%
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [chipod, cham] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
        
        eps_cham_all = [eps_cham_all(:) ; cham.eps(:) ] ;
        chi_cham_all = [chi_cham_all(:) ; cham.chi(:) ] ;
        N2_cham_all  = [N2_cham_all(:)  ; cham.N2(:) ] ;
        Tz_cham_all  = [Tz_cham_all(:)  ; cham.Tz(:) ] ;
        
        eps_chi_all = [eps_chi_all(:) ; chipod.eps(:) ] ;
        chi_chi_all = [chi_chi_all(:) ; chipod.chi(:) ] ;
        N2_chi_all  = [N2_chi_all(:)  ; chipod.N2(:) ] ;
        Tz_chi_all  = [Tz_chi_all(:)  ; chipod.Tz(:) ] ;
        
    end % idx
    
    subplot(211)
    histogram(  log10(chi_chi_all./chi_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    xlim([-2 2])
    %ylim([0 1.5])
    grid on
    hold on
    freqline(nanmean(log10(chi_chi_all./chi_cham_all)),cols(iax))
    text(1.25,0.8-(iax*0.1),['\mu= ' num2str(roundx(nanmean(log10(chi_chi_all./chi_cham_all)),2))],'color',cols(iax),'fontsize',14)
    hold on
    
    subplot(212)
    hh = histogram(  log10(eps_chi_all./eps_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    xlim([-3 3])
    ylim([0 1.2])
    grid on
    hold on
    freqline(nanmean(log10(eps_chi_all./eps_cham_all)),cols(iax))
    text(1.5,0.8-(iax*0.1),['\mu= ' num2str(roundx(nanmean(log10(eps_chi_all./eps_cham_all)),2))],'color',cols(iax),'fontsize',14)
    hold on
    
    h   = [h hh];
    iax = iax+1
    
end % dp
%
subplot(212)
legend(h,'1 profile','10 profiles','50 profiles')
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]','fontsize',16)
ylabel('pdf','fontsize',16)

subplot(211)
xlabel('log_{10}[\chi_{\chi}/\chi]','fontsize',16)
ylabel('pdf','fontsize',16)

print( fullfile(fig_dir,[project_short '_eps_ratio_hist_diff_prof_avg_' 'Pmin_' num2str(Pmin) MakeChiPathStr(Params)]), '-dpng')

%%