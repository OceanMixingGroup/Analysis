%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_Overview_Plots.m
%
% Making plots for  eq14 chi-pod overview notes
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

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

Params.gamma     = 0.2;
Params.fmax      = 7  ;
Params.z_smooth  = 1  ;
Params.resp_corr = 0  ;
Params.fc        = 99 ;

dz = 2 ;
cnums_to_get = get_cham_cnums_eq14;
cnums_to_get = cnums_to_get(cnums_to_get>500);
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

screen_chi = 1 ;
screen_ml  = 0 ;
Pmin       = 0 ;

eq14_patches_paths


%%
%[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
if exist(fullfile(analysis_dir,project_short,'data',save_name),'file')==2
    load(fullfile(analysis_dir,project_short,'data',save_name))
else
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat')


%% Pcolor of chipod & cham chi, and N2,Tz

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ; cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hold on
plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
plot(zml_cnum,zml,'k')
caxis([-6 -2])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
plot(zml_cnum,zml,'k')
caxis([-4 -0])
colorbar
ylabel('P [db]','fontsize',16)
xlabel('cast #','fontsize',16)
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])%%
figname = [project_short '_Pcolor_BothChi_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir, figname), '-dpng')


%% Pcolor of chipod & cham *eps*, and N2,Tz

figure(2);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1) ;
ezpc(cham.cnum,cham.P,log10(cham.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chameleon')
ylabel('P [db]','fontsize',16)

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chi-pod')
ylabel('P [db]','fontsize',16)

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
plot(zml_cnum,zml,'k')
caxis([-6 -2])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
plot(zml_cnum,zml,'k')
caxis([-4 -0])
colorbar
ylabel('P [db]','fontsize',16)
xlabel('cast #','fontsize',16)
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])
%
figname = [project_short '_Pcolor_BothEps_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir, figname), '-dpng')



%% Combine above 2 plots - plot chi,eps,n2,dtdz on 6X1 panel plot

figure(1);clf
agutwocolumn(1)
wysiwyg

rr=6 ;
cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3) ;
ezpc(cham.cnum,cham.P,log10(cham.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chameleon')
ylabel('P [db]','fontsize',16)

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,log10(chipod.eps))
hold on
plot(zml_cnum,zml,'k')
caxis([-11 -4])
colorbar
title('log_{10} \epsilon chi-pod')
ylabel('P [db]','fontsize',16)

ax5 = subplot(rr,cc,5);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
plot(zml_cnum,zml,'k')
caxis([-6 -2])
colorbar
ylabel('P [db]','fontsize',16)
title('log_{10} N^2')

ax6 = subplot(rr,cc,6);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
plot(zml_cnum,zml,'k')
caxis([-4 -0])
colorbar
ylabel('P [db]','fontsize',16)
xlabel('cast #','fontsize',16)
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6])

figname = [project_short '_Pcolor_Both_epsANDChi_N2_Tz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir, figname), '-dpng')


%% Reload data, screening surface and mixed layer regions

clear cham chipod

dz = 2 ;
screen_chi = 1
Pmin       = 20
screen_ml  = 1

save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
if exist(fullfile(analysis_dir,project_short,'data',save_name),'file')==2
    load(fullfile(analysis_dir,project_short,'data',save_name))
else
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end


%% Plot ratio of chipod/cham chi and eps, vs. actual epsilon

figure(1);clf
h = figure ; clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(211) ;
histogram2( log10(cham.eps),log10(chipod.chi ./ cham.chi ), 'DisplayStyle','tile','Normalization','pdf')
hold on
ylim([-2.5 2.5])
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\chi_{\chi}/\chi]','fontsize',16)
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])
hline(0,'k')
caxis([0 0.5])

ax2 = subplot(212);
histogram2( log10(cham.eps), log10(chipod.eps ./ cham.eps),50, 'DisplayStyle','tile','Normalization','pdf')
hold on
xlim([-8.5 -4.5])
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\epsilon_{\chi}/\epsilon]','fontsize',16)
hline(0,'k')
caxis([0 0.4])

figname = [project_short '_ratios_vs_eps_' num2str(screen_chi) '_Pmin_' num2str(Pmin) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir, figname), '-dpng')

%% Plot CHAMELEON eps vs chi, normalized so slope is equal to 1/2*gamma

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
title([project_short ' Chameleon ' num2str(dz) 'm binned'])

%fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_Pmin_' num2str(Pmin)]
fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_' MakeChiPathStr(Params)]
print( fullfile(fig_dir,fname),'-dpng')


%% 2D histograms of chi, eps vs chameleon

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

% Histograms of ratio of chi-pod epsilon to chameleon epsilon

figure(13);clf
agutwocolumn(0.8)
wysiwyg

subplot(211)

h1 = histogram(log10( chipod.chi(:) ./ cham.chi(:) ),'EdgeColor','none','Normalization','pdf') ;
hold on
xlim([-2 2])
grid on
xlabel('log_{10}[\chi_{\chi} /\chi ]','fontsize',16)
ylabel('pdf','fontsize',16)
%title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])
%title([project_short ' ' num2str(dz) ' m binned '])
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])
freqline( nanmean(log10( chipod.chi(:) ./ cham.chi(:) )))
text(1,0.4,['\mu= ' num2str( roundx(nanmean(log10( chipod.chi(:) ./ cham.chi(:) )),2))],'fontsize',15)

subplot(212)
h2 = histogram(log10( chipod.eps(:) ./ cham.eps(:) ),'EdgeColor','none','Normalization','pdf') ;
hold on
xlim([-3 3])
grid on
xlabel('log_{10}[\epsilon_{\chi} /\epsilon ]','fontsize',16)
ylabel('pdf','fontsize',16)
%title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])
%title([project_short ' ' num2str(dz) ' m binned '])
freqline( nanmean(log10( chipod.eps(:) ./ cham.eps(:) )))
text(1,0.4,['\mu= ' num2str( roundx(nanmean(log10( chipod.eps(:) ./ cham.eps(:) )),2))],'fontsize',15)


clear figname
figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_screen_chi_' num2str(screen_chi) '_screenml_' num2str(screen_ml)  '_' MakeChiPathStr(Params)]
print( fullfile(fig_dir,figname),'-dpng')
clear figname

%% 2D hist of the ratio of chipod/cham vs depth?

%
P = repmat(chipod.P,1,length(chipod.cnum));

figure(33);clf
agutwocolumn(1)
wysiwyg

subplot(211)
histogram2( real(log10( chipod.chi ./ cham.chi )), P, 'DisplayStyle','tile')
axis ij
ylabel('P[db]','fontsize',16)
xlabel('log_{10}[\chi_{\chi}/\chi]','fontsize',16)
xlim([-3 3])
colorbar
caxis([0 1000])
freqline(0,'k-')
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ', dz=' num2str(dz)])


subplot(212)
histogram2( real(log10( chipod.eps ./ cham.eps )), P, 'DisplayStyle','tile','Edgecolor','none')
axis ij
ylabel('P[db]','fontsize',16)
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]','fontsize',16)
xlim([-3 3])
caxis([0 1000])
colorbar
freqline(0,'k-')

figname=['eq14_chi_eps_Vs_P_2Dhist_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')


%% plot chi vs chi and eps vs eps for different depth bin averaging

clear cham chipod
%close all

%addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

screen_chi = 1 ;
screen_ml  = 1 ;
Pmin       = 20;

zmin=0  ;
zmax=200;

% addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
% addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/
%
% eq14_patches_paths

cnums_to_get = get_cham_cnums_eq14 ;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);

figure(8);clf
agutwocolumn(1)
wysiwyg

iax=1 ; rr=3 ; cc=2 ;
for dz = [2 10 50]
    
    clear chipod cham
    %    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,Pmin,screen_chi,screen_ml) ;
    
    save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
    if exist(fullfile(analysis_dir,project_short,'data',save_name),'file')==2
        load(fullfile(analysis_dir,project_short,'data',save_name))
    else
        [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
    end
    
    subplot(rr,cc,iax)
    hh = histogram2(  real(log10(cham.chi)),log10(chipod.chi),80,'DisplayStyle','tile','XBinEdges',[-11:0.2:-4],'YBinEdges',[-11:0.2:-4],'Edgecolor','none');
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
    hh = histogram2(  real(log10(cham.eps)),log10(chipod.eps),80,'DisplayStyle','tile','XBinEdges',[-9:0.2:-4],'YBinEdges',[-9:0.2:-4],'Edgecolor','none');
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
figname = ['eq14_chiVscham_chiANDeps_diff_dz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')
clear figname

%% histogram of epsilon ratio for different size depth averaging

clear cham chipod

screen_chi = 1 ;
Pmin       = 20;
screen_ml  = 1 ;

zmin=0  ;
zmax=200;

cnums_to_get = get_cham_cnums_eq14 ;
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);


figure(9);clf
agutwocolumn(1)
wysiwyg

iax=1
cols=['b','r','g']
h=[]

for dz = [2 10 50]
    
    clear chipod cham
    %    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,Pmin,screen_chi,screen_ml) ;
    save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
    if exist(fullfile(analysis_dir,project_short,'data',save_name),'file')==2
        load(fullfile(analysis_dir,project_short,'data',save_name))
    else
        [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
    end
    
    subplot(2,1,1)
    hh=histogram( log10( chipod.chi(:) ./ cham.chi(:)),[-3:0.1:3], 'Normalization','pdf','DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax));
    grid on
    xlim([-2 2])
    ylim([0 1.4])
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
    epsrat=nanmean(log10( chipod.eps(:) ./ cham.eps(:)));
    text(1,0.8-(iax*0.1),['\mu= ' num2str(roundx(epsrat,2))],'fontsize',15,'color',cols(iax))
    
    %h=[h hh];
    
    iax=iax+1;
    
end


xlabel(['log_{10}[\epsilon_{\chi}/\epsilon]'],'fontsize',16)
ylabel('pdf')
%
subplot(211)
legend(h,'2m','10m','50m')
xlabel(['log_{10}[\chi_{\chi}/\chi]'],'fontsize',16)
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ])
ylabel('pdf')

%
clear figname
figname=['eq14_chiVscham_hist_diff_dz_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
print(fullfile(fig_dir,figname),'-dpng')
clear figname


%% Plot chi vs chi, eps vs eps, for different # profiles averaged

clear cham chipod

screen_chi= 1 ;
screen_ml = 1 ;
Pmin      = 20;

dz        = 10; % bin size

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(11);clf
agutwocolumn(1)
wysiwyg

iax=1;
for dp = [2 10 20]
    
    eps_cham_all = [];
    chi_cham_all = [];
    N2_cham_all = [];
    Tz_cham_all = [];
    
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
    hh = histogram2(  real(log10(chi_cham_all)),log10(chi_chi_all),40,'DisplayStyle','tile','XBinEdges',[-11:0.2:-4],'YBinEdges',[-11:0.2:-4],'Edgecolor','none');
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
    hh = histogram2(  real(log10(eps_cham_all)),log10(eps_chi_all),20,'DisplayStyle','tile','XBinEdges',[-8.5:0.2:-4],'YBinEdges',[-8.5:0.2:-4],'Edgecolor','none');
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
figname=['eq14_chiVscham_chiANDeps_diff_prof_avg_screen_chi_' num2str(screen_chi) '_' MakeChiPathStr(Params)]
%figname=['eq14_chiVscham_chiANDeps_diff_prof_avg_screen_chi_' num2str(screen_chi) '_Pmin_80']
print(fullfile(fig_dir,figname),'-dpng')


%% May 1 2017 - plot histograms of eps_chi/eps for differnt # prof avg.

clear cham chipod
% close all

% Params.gamma    = 0.2;
% Params.fmax     = 7  ;
% Params.z_smooth = 1  ;
% Params.resp_corr= 0  ;    % correct TP spectra for freq response of thermistor?
% Params.fc       = 99 ;    % cutoff frequency for response correction

screen_chi = 1 ;
Pmin       = 20;
screen_ml  = 1 ;

dz = 10 % bin size

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(12);clf
agutwocolumn(1)
wysiwyg

h=[] ; iax=1 ;cols=['b','r','g'];
for dp = [2 10 20 ]
    
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
    hh = histogram(  log10(chi_chi_all./chi_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax)) ;
    xlim([-2 2])
    ylim([0 1.5])
    grid on
    hold on
    freqline(nanmean(log10(chi_chi_all./chi_cham_all)),cols(iax))
    text(1,0.8-(iax*0.1),['\mu= ' num2str(roundx(nanmean(log10(chi_chi_all./chi_cham_all)),2))],'color',cols(iax),'fontsize',14)
    hold on
    
    subplot(212)
    histogram(  log10(eps_chi_all./eps_cham_all),[-2:0.15:2],'Normalization','pdf','FaceAlpha',0.5,'DisplayStyle','stair','LineWidth',2,'EdgeColor',cols(iax)) ;
    xlim([-3 3])
    ylim([0 1.2])
    grid on
    hold on
    freqline(nanmean(log10(eps_chi_all./eps_cham_all)),cols(iax))
    text(1.5,0.8-(iax*0.1),['\mu= ' num2str(roundx(nanmean(log10(eps_chi_all./eps_cham_all)),2))],'color',cols(iax),'fontsize',14)
    hold on
    
    h   = [h hh];
    iax = iax+1 ;
    
end % dp
%
subplot(212)
xlabel('\epsilon_{\chi}/\epsilon','fontsize',16)
ylabel('pdf','fontsize',16)

subplot(211)
legend(h,'2 profile','10 profiles','20 profiles')
xlabel('\chi_{\chi}/\chi','fontsize',16)
ylabel('pdf','fontsize',16)
title(['fmax=' num2str(Params.fmax) ', zsmooth= ' num2str(Params.z_smooth) ])

print( fullfile(fig_dir,['eq14_eps_ratio_hist_diff_prof_avg_' MakeChiPathStr(Params)]), '-dpng')


%
%% See if gamma computed from multi-profile averages of N2,Tz,chi,eps is 0.2?
% compare to gamma computed from individual 1m data points in every profile

clear ; %close all

dz = 10; % bin size to average over

eq14_patches_paths

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

cnums = get_cham_cnums_eq14;

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
title(['1mavg'])%, profiles ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])

ax2 = subplot(122) ;
boxplot(log10(gam_avg))
grid on
hline(log10(0.2),'k--')
title(['profile-averaged, ' num2str(dz) ' m binned'])

linkaxes([ax1 ax2])

figname=[project_short '_gamma_point_avg_box_' num2str(dz) 'mbinned']
print(fullfile(fig_dir,figname),'-dpng')

%%