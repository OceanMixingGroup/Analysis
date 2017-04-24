%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_eps_profiles.m
%
% Plot profiles of epsilon from chameleon and chi-pod method for EQ14.
% Compare time-average profiles for groups of casts?
%
%
%----------------
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot single profiles: chameleon, chi-pod binned and patch

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;

dz=10

eq14_patches_paths
figdir1 = fullfile( fig_dir, 'chi_profiles', ['fmax_' num2str(Params.fmax) '_zsmooth_' num2str(Params.z_smooth)]);
ChkMkDir(figdir1)
figdir2 = fullfile( fig_dir, 'eps_profiles', ['fmax_' num2str(Params.fmax) '_zsmooth_' num2str(Params.z_smooth)]);
ChkMkDir(figdir2)

for cnum=1:100:3000
    try
        h = PlotChiProfileCompare_eq14(cnum,Params,dz)
        print( fullfile( figdir1, ['eq14_profile_' num2str(cnum) '_chi_profiiles_compare'] ),'-dpng')

        h = PlotEpsProfileCompare_eq14(cnum,Params,dz)
        print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_eps_profiiles_compare'] ),'-dpng')
        pause(1)
    end
end

%% do same, but plot average of multiple profiles

clear ; close all

Params.gamma = 0.05;
Params.fmax  = 7  ;
Params.z_smooth=10;

dz=10

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

dp=40
Pmin=0;

figdir2 = fullfile( fig_dir, ['chi_eps_profiles_' num2str(dp) 'profavgs'],['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128']);
ChkMkDir(figdir2)


for cnum=1:50:3000
    try

        % avg +/- dp profiles
        cnums_to_get = (cnum-dp/2) : (cnum+dp/2);   
        
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin);

        figure(1);clf
        agutwocolumn(1)
        wysiwyg
                       
        clear chb avg
        % regular chi-pod method on binned data
        load( fullfile(path_chipod_bin,['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
        chb=avg;clear avg
        
        % chamelon data
        load(fullfile(path_cham_avg,['EQ14_' sprintf('%04d',cnum) '_avg.mat']))

        subplot(221)
        plot(log10(chb.chi1),chb.P,'color',0.6*[1 1 1])
        hold on
        plot(log10(avg.CHI),avg.P,'k')
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        xlabel('log_{10}[\chi]')
        title(['profile ' num2str(cnum)])
        
        
        subplot(222)
        hcham=plot(log10(chi_cham_avg),P_cham,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chi_chi_avg),P_chi,'d-','color',0.6*[1 1 1],'linewidth',2);
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        legend([hcham hchi],'cham','\chi pod','location','best')
        xlabel('log_{10}[\chi]')
        title(['profiles ' num2str(cnum-dp) ' - ' num2str(cnum+dp)])

        subplot(223)
        plot(log10(chb.eps1),chb.P,'color',0.6*[1 1 1])
        hold on
        plot(log10(avg.EPSILON),avg.P,'k')
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        xlabel('log_{10}[\epsilon]')
        title(['profile ' num2str(cnum)])
        
        
        subplot(224)
        hcham=plot(log10(eps_cham_avg),P_cham,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(eps_chi_avg),P_chi,'d-','color',0.6*[1 1 1],'linewidth',2);
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        legend([hcham hchi],'cham','\chi pod','location','best')
        xlabel('log_{10}[\epsilon]')
        title(['profiles ' num2str(cnum-dp) ' - ' num2str(cnum+dp)])
        
        print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_eps_profiiles_compare'] ),'-dpng')
        
        pause(0.1)
    %catch
    end
end



%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz = 10 % bin size
zmin=0;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths

cnums_to_get = get_cham_cnums_eq14 ;

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax)


% Nan values in mixed layer

Pmin = 80;

clear ib
ib = find(cham.P<Pmin);
cham.eps(ib) = nan;

clear ib
ib = find(chipod.P<Pmin);
chipod.eps(ib) = nan;

%
%% Histograms of ratio of chi-pod epsilon to chameleon epsilon (10bins)

figure(1);clf
h1 = histogram(log10( chipod.eps(:) ./ cham.eps(:) ),'EdgeColor','none','Normalization','pdf');
hold on
%h2 = histogram(log10(patch./cham),'EdgeColor','none','Normalization','pdf');
xlim([-4 4])
grid on
%legend([h1 h2],'bin/cham','patch/cham')
xlabel('log_{10}[\epsilon_{\chi} /\epsilon ]')
ylabel('pdf')
title([project_short ' ' num2str(dz) ' m binned, Pmin=' num2str(Pmin)])

%

figname=[project_short '_' num2str(dz) 'mbinned_eps_ratios_Pmin' num2str(Pmin)]
print( fullfile(fig_dir,figname),'-dpng')
SetNotesFigDir
print( fullfile(NotesFigDir,figname),'-dpng')


% Plot eps vs chi, normalized so slope is equal to 1/2*gamma

figure(2);clf
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


%% Now want to do same, but average some profiles together to see if
% converges to gamma=0.2

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth =10 ;

dz = 10 % bin size
Pmin = 80

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg

% figure(2);clf
% agutwocolumn(1)
% wysiwyg

iax=1;
for dp=[10 50 100 500]
    
    normx_all = [];
    normy_all = [];
    
    normx_all_chi = [];
    normy_all_chi = [];
    
    for ix = 1:round(3000/dp)%30
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        
        %clear eps_cham chi_cham N2_cham Tz_cham
        %clear eps_chi chi_chi N2_chi Tz_chi
        %    [eps_cham, chi_cham, N2_cham, Tz_cham, eps_chi, chi_chi, N2_chi, Tz_chi] =...
        %        Get_binned_data(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short);
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_cham,P_chi] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin);
        
        normy_all = [ normy_all(:) ; chi_cham_avg./(Tz_cham_avg.^2)] ;
        normx_all = [ normx_all(:) ; eps_cham_avg./N2_cham_avg ];
        
        normy_all_chi = [ normy_all_chi(:) ; chi_chi_avg./(Tz_chi_avg.^2)] ;
        normx_all_chi = [ normx_all_chi(:) ; eps_chi_avg./N2_chi_avg ];
        
    end % idx
    
    figure(1)
    subplot(2,2, iax)
    %hh=histogram2( log10(normx_all) , real(log10(normy_all)),80,'DisplayStyle','tile')
    %hh=histogram2( log10(chi_chi./(Tz_chi.^2)) , real(log10(eps_chi./N2_chi)),80,'DisplayStyle','tile')
    hh=scatter( log10(normx_all) , log10(normy_all),'filled','MarkerFaceAlpha',0.1)
    %loglog(chi_cham./(Tz_cham.^2),eps_cham./N2_cham,'.','color',0.75*[1 1 1])
    grid on
    hold on
    xvec=linspace(1e-7,1e0,100);
    h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
    h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
    h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
    ylim([-7.5 -1])
    xlim([-5.5 -1])
    ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
    xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
    legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
    title([project_short ' 10m binned ,' num2str(dp) ' profile averages'])
    
    %     fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_' num2str(dp) 'profileAvg']
    %     print( fullfile(fig_dir,fname),'-dpng')
    
    %         figure(2)
    %     subplot(2,2, iax)
    %     %hh=histogram2( log10(normx_all) , real(log10(normy_all)),80,'DisplayStyle','tile')
    %     %hh=histogram2( log10(chi_chi./(Tz_chi.^2)) , real(log10(eps_chi./N2_chi)),80,'DisplayStyle','tile')
    %     hh=scatter( log10(normx_all_chi) , log10(normy_all_chi),'filled','MarkerFaceAlpha',0.1)
    %     %loglog(chi_cham./(Tz_cham.^2),eps_cham./N2_cham,'.','color',0.75*[1 1 1])
    %     grid on
    %     hold on
    %     xvec=linspace(1e-7,1e0,100);
    %     h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
    %     h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
    %     h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
    %     ylim([-7.5 -1])
    %     xlim([-5.5 -1])
    %     ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
    %     xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
    %     legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
    %     title([project_short ' 10m binned ,' num2str(dp) ' profile averages'])
    %
    
    iax=iax+1;
    
end % dp (# profiles averaged together)
%
figure(1);shg
fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_diffNprofileAvg_Pmin' num2str(Pmin)]
print( fullfile(fig_dir,fname),'-dpng')

%% plot locations of chameleon casts

clear ; close all

% first find profiles near each other
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

%
figure(1);clf
subplot(211)
plot(cham.castnumber,cham.lon)
subplot(212)
plot(cham.castnumber,cham.lat)

%% plot groups of profiles averaged together in depth bins

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg

for whcase=1:6
    
    switch whcase
        case 1
            cnum_range = [0 500];
        case 2
            cnum_range = [500 1000];
        case 3
            cnum_range = [1000 1500];
        case 4
            cnum_range = [1500 2000];
        case 5
            cnum_range = [2000 2500];
        case 6
            cnum_range = [2500 3000];
    end
    
    clear cnums
    cnums_to_get = [cnum_range(1) : cnum_range(2) ];
    
    Pmin=0;
    
    [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
        Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin)
        
    screen=1
    if screen==1
        % try screening out some spikes in chipod data that give huge epsilons
        clear ib
        ib=find( medfilt1(Tz_chi_avg,5) ./ Tz_chi_avg  >2 ) ;
        eps_chi_avg(ib)=nan;
        
        clear ib
        ib = find(log10(N2_chi_avg)>-2.5);
        eps_chi_avg(ib)=nan;
                
        eps_chi_avg(find(log10(eps_chi_avg)>-4))=nan;
        eps_chi_avg(find(log10(eps_chi_avg)<-8.5))=nan;
    end    
    
    %
    figure(1)
    subplot(2,3,whcase)
    hcham = plot(log10(eps_cham_avg),P_cham,'ms-','linewidth',2) ;
    hold on
    hchi = plot(log10(eps_chi_avg),P_chi,'rp-','linewidth',2) ;
    axis ij
    grid on
    xlim([-9.5 -4])
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    if whcase==5
        legend([ hcham hchi],'cham','chipod','location','best')
    end
    title([num2str(cnum_range(1)) '-' num2str(cnum_range(2)) ])
    
end

%%

figure(1)
eq14_patches_paths
figname=[project_short '_eps_prof_comparisons_' num2str(dz) 'mbins']
print( fullfile(fig_dir,figname),'-dpng')


%% same as above but different method

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg


cnum_range = [2600 2620];

clear cnums
cnums_to_get = [cnum_range(1) : cnum_range(2) ];

Pmin=0;

[chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
    path_cham_avg,dz,Params,cnums_to_get,project_short,0,200);

%

[eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
[eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);

figure(1);clf
agutwocolumn(0.75)
wysiwyg

subplot(121)
plot(log10(chipod.eps),chipod.P,'.','color',0.8*[1 1 1])
hold on
plot(log10(cham.eps),cham.P,'.','color',0.2*[1 1 1])
hcham=plot(log10(eps_cham_avg),p_cham,'linewidth',2,'color',0.2*[1 1 1])
hold on
hchi=plot(log10(eps_chi_avg),p_chi,'linewidth',2,'color',0.8*[1 1 1])
axis ij
grid on
legend([hcham hchi],'cham','\chi pod','location','best')
ylim([0 200])
xlim([-12 -2])
title(['profiles ' num2str(cnum_range(1)) ' - ' num2str(cnum_range(2))])
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('P','fontsize',16)

% Try screening spikes in chipod data and re-plot?

% try screening out some spikes in chipod data that give huge epsilons
%     clear ib
%     ib=find( medfilt1(Tz_chi,5) ./ Tz_chi  >2 ) ;
%     eps_chi(ib)=nan;

%     clear ib
%     ib = find(log10(N2_chi)>-2.5);
%     eps_chi(ib)=nan;

%     clear ib
%     ib = find(log10(cham.EPSILON)<-8.5);
%     cham.EPSILON(ib) = nan;

chipod.eps(find(log10(chipod.eps)>-4))=nan;

chipod.eps(find(log10(chipod.eps)<-8.5))=nan;

[eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
[eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);

% figure(2);clf
% agutwocolumn(0.75)
% wysiwyg
subplot(122)
plot(log10(chipod.eps),chipod.P,'.','color',0.8*[1 1 1])
hold on
plot(log10(cham.eps),cham.P,'.','color',0.2*[1 1 1])
hcham=plot(log10(eps_cham_avg),p_cham,'linewidth',2,'color',0.2*[1 1 1])
hold on
hchi=plot(log10(eps_chi_avg),p_chi,'linewidth',2,'color',0.8*[1 1 1])
axis ij
grid on
legend([hcham hchi],'cham','\chi pod','location','best')
ylim([0 200])
xlim([-12 -2])
%title(['profiles ' num2str(cnum_range(1)) ' - ' num2str(cnum_range(2))])
title('after screening')
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('P','fontsize',16)

figname = [project_short '_profiles_' num2str(cnum_range(1)) '_' num2str(cnum_range(2)) '_eps_profiles_avg']
print(fullfile(fig_dir,figname), '-dpng')

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

% [N2_bin  zout Nobs] = binprofile(cham.N2(:,iCham)           ,cham.P(:,iCham), 0, dz, 200,1);
% [Tz_bin  zout Nobs] = binprofile(cham.DTDZ_RHOORDER(:,iCham),cham.P(:,iCham), 0, dz, 200,1);
% [chi_bin zout Nobs] = binprofile(cham.CHI(:,iCham)          ,cham.P(:,iCham), 0, dz, 200,1);
% [eps_bin zout Nobs] = binprofile(cham.EPSILON(:,iCham)      ,cham.P(:,iCham), 0, dz, 200,1);

[eps_bin zout Nobs] = binprofile(eps ,P, 0, dz, 200,1);
[chi_bin zout Nobs] = binprofile(chi ,P, 0, dz, 200,1);
[N2_bin  zout Nobs] = binprofile(N2  ,P, 0, dz, 200,1);
[Tz_bin  zout Nobs] = binprofile(Tz  ,P, 0, dz, 200,1);

gam_avg = ComputeGamma(N2_bin,Tz_bin,chi_bin,eps_bin);

figure(1);clf
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


%% Try varying # of profiles averaged to see how many it takes to converge

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2 ;
Params.fmax = 7 ;

% patch parameters
patch_size_min = 0.4 ;
usetemp = 1 ;
minR2   = 0.0 ;

eq14_patches_paths

dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches');

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])

figure(1);clf
agutwocolumn(1)
wysiwyg

e1_all = [] ;
e2_all = [] ;
cham_all = [] ;

prof_start=2200
prof_extra = [5 25 50 100 150 300 ];
for whcase=1:6
    
    switch whcase
        case 1
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
        case 2
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
        case 3
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
        case 4
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
        case 5
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
        case 6
            cnum_range = [prof_start prof_start+prof_extra(whcase)];
    end
    
    clear cnums
    cnums = [cnum_range(1) : cnum_range(2) ];
    
    e1=[]; P1=[];
    e2=[]; P2=[]; N2=[] ; Tz=[] ; tpvar=[]; fspd=[]; cnall=[];
    
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            
            % patch chi-pod
            clear avg
            load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
            e1 = [e1(:) ; avg.eps1(:)];
            P1 = [P1(:) ; avg.P(:)   ];
            
            % binned chi-pod
            clear avg
            load(fullfile(path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
            e2 = [e2(:) ; avg.eps1(:)];
            P2 = [P2(:) ; avg.P(:)   ];
            N2 = [N2(:) ; avg.N2(:) ] ;
            Tz = [Tz(:) ; avg.dTdz(:) ];
            tpvar = [tpvar(:) ; avg.TP1var(:) ];
            fspd = [fspd(:) ; avg.fspd(:) ];
            cnall = [ cnall(:) ; repmat(cnum,length(avg.P),1)];
            
        end
    end
    
    
    iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));
    
    e2(find(log10(e2)>-4))=nan;
    
    [dataout1 zout1 Nobs] = binprofile(e1,P1, 0, 10, 200,1);
    [dataout2 zout2 Nobs] = binprofile(e2,P2, 0, 10, 200,1);
    [cham_bin zout_cham Nobs] = binprofile(cham.EPSILON(:,iCham),cham.P(:,iCham), 0, 10, 200,1);
    
    e2_all = [e2_all dataout2] ;
    cham_all = [cham_all cham_bin ];
    
end

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

for iax=1:6
    
    subplot(2,3,iax)
    plot(log10(e2_all(:,iax)),zout2,'ro-')
    hold on
    plot(log10(cham_all(:,iax)),zout_cham,'k')
    axis ij
    grid on
    xlim([-11 -4])
    title([num2str(prof_start) ' - ' num2str(prof_start+prof_extra(iax))])
end

print(fullfile(fig_dir,[project_short '_eps_prof_diffN']),'-dpng')

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

for iax=1:6
    
    subplot(2,3,iax)
    histogram(log10(e2_all(:,iax)./cham_all(:,iax)),[-2.5:0.25:2.5],'Normalization','pdf')
    hold on
    grid on
    xlim(2.5*[-1 1])
    ylim([0 1])
    title([num2str(prof_extra(iax))])
    nanstd(e2_all(:,iax)./cham_all(:,iax))
    
end

%%