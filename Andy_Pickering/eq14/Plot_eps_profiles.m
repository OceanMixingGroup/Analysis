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

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;
screen_chi=0

dz=10

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

dp = 10
Pmin = 20 ;
screen_ml=1

figdir2 = fullfile( fig_dir, ['chi_eps_profiles_' num2str(dp) 'profavgs'],['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_screen_chi_' num2str(screen_chi)]);
ChkMkDir(figdir2)

for cnum=551:50:3000
    try

        % avg +/- dp profiles
        cnums_to_get = (cnum-dp/2) : (cnum+dp/2);  
        
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get = setdiff(cnums_to_get,bad_prof);
        
        [chipod,cham] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);

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
        hline(80,'k--')
        
        subplot(222)
        hcham=plot(log10(cham.chi),cham.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod.chi),chipod.P,'d-','color',0.6*[1 1 1],'linewidth',2);
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        legend([hcham hchi],'cham','\chi pod','location','best')
        xlabel('log_{10}[\chi]')
        title(['profiles ' num2str(cnum-dp) ' - ' num2str(cnum+dp)])
        hline(80,'k--')
        
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
                hline(80,'k--')
        
        subplot(224)
        hcham=plot(log10(cham.eps),cham.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod.eps),chipod.P,'d-','color',0.6*[1 1 1],'linewidth',2);
        axis ij
        grid on
        xlim([-12 -3])
        ylim([0 200])
        legend([hcham hchi],'cham','\chi pod','location','best')
        xlabel('log_{10}[\epsilon]')
        title(['profiles ' num2str(cnum-dp) ' - ' num2str(cnum+dp)])
                hline(80,'k--')
                
        %print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_eps_profiiles_compare'] ),'-dpng')
        
        pause(0.1)
    %catch
    end
end


%% Make similar profile plots, but plot average of different # profiles on same figure?

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;


dz=10

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq14_patches_paths

dp = 100
Pmin = 20 ;
screen_ml=1
screen_chi=1

figdir2 = fullfile( fig_dir, ['chi_eps_profiles_diffNprof_profavgs'],['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_screen_chi_' num2str(screen_chi)]);
ChkMkDir(figdir2)

rr=2; cc=3;

for cnum=551:25:3000
    try

        figure(1);clf
        agutwocolumn(1)
        wysiwyg
        
        clear chb avg
        % regular chi-pod method on binned data
        load( fullfile(path_chipod_bin,['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
        chb=avg;clear avg
        chb.eps1(find(log10(chb.eps1)<-8.5))=nan;
        chb.P(find(chb.P<Pmin))=nan;
        chb = discard_convection_eq14_chi(chb);
        
        % chamelon data
        load(fullfile(path_cham_avg,['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
        avg.EPSILON(find(log10(avg.EPSILON)<-8.5))=nan;
        avg = discard_convection_eq14_cham(avg);
        avg.P(find(avg.P<Pmin))=nan;
        
        ax1=subplot(rr,cc,1);
        hcham=plot(log10(avg.CHI),avg.P,'k.-')
        hold on
        hchi=plot(log10(chb.chi1),chb.P,'b.-')%'color',0.6*[1 1 1])        
        axis ij
        grid on
        xlim([-12 -4])
        ylim([0 200])
        xlabel('log_{10}[\chi]')
        %legend([hcham hchi],'cham','\chi pod','location','northwest','orientation','horizontal')
        title(['profile ' num2str(cnum)])
        %hline(80,'k--')
        
        ax2=subplot(rr,cc,cc+1);
        plot(log10(avg.EPSILON),avg.P,'k.-')
        hold on
        plot(log10(chb.eps1),chb.P,'b.-')
        hold on       
        axis ij
        grid on
        xlim([-9 -4])
        ylim([0 200])
        xlabel('log_{10}[\epsilon]')
        title(['profile ' num2str(cnum)])

        
        dp = 10;
        % avg +/- dp profiles
        cnums_to_get_1 = (cnum-dp/2) : (cnum+dp/2);  
        
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get_1 = setdiff(cnums_to_get_1,bad_prof);
        
        [chipod , cham] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get_1,project_short,Pmin,screen_chi,screen_ml);

        dp = 50;
        % avg +/- dp profiles
        cnums_to_get_2 = (cnum-dp/2) : (cnum+dp/2);  
        
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get_2 = setdiff(cnums_to_get_2,bad_prof);
        
        [chipod2 , cham2] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get_2,project_short,Pmin,screen_chi,screen_ml);
        
        subplot(rr,cc,2)
        hcham=plot(log10(cham.chi),cham.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod.chi),chipod.P,'b-','linewidth',1);
        axis ij
        grid on
        xlim([-10 -4])
        ylim([0 200])
        xlabel('log_{10}[\chi]')
        title(['profiles ' num2str(cnums_to_get_1(1)) ' - ' num2str(cnums_to_get_1(end))])

        subplot(rr,cc,3)
        hcham=plot(log10(cham2.chi),cham2.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod2.chi),chipod2.P,'b-','linewidth',1);
        axis ij
        grid on
        xlim([-10 -4])
        ylim([0 200])
        xlabel('log_{10}[\chi]')
        title(['profiles ' num2str(cnums_to_get_2(1)) ' - ' num2str(cnums_to_get_2(end))])

        subplot(rr,cc,5)
        hcham=plot(log10(cham.eps),cham.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod.eps),chipod.P,'b-','linewidth',2);
        axis ij
        grid on
        xlim([-8.5 -4])
        ylim([0 200])
        xlabel('log_{10}[\epsilon]')
        title('cham = black')
                
        subplot(rr,cc,6)
        hcham=plot(log10(cham2.eps),cham2.P,'ko-','linewidth',2);
        hold on
        hchi = plot(log10(chipod2.eps),chipod2.P,'b-','linewidth',2);
        axis ij
        grid on
        xlim([-8.5 -4])
        ylim([0 200])
        xlabel('log_{10}[\epsilon]')
        title('\chi pod = blue','color','b')
        
        print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_eps_profiles_compare'] ),'-dpng')
        
       % pause(0.1)
    catch
    end
end % cnum





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

screen_chi = 1

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
    cnums_to_get = [cnum_range(1) : cnum_range(2) ] ;

    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);

    Pmin=0;
    
    [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
        Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi)
        
    screen=0
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

%

figure(1)
eq14_patches_paths
%figname=[project_short '_eps_prof_comparisons_' num2str(dz) 'mbins']
figname=[project_short '_eps_prof_comparisons_' num2str(dz) 'mbins_screen_chi_' num2str(screen_chi)]
print( fullfile(fig_dir,figname),'-dpng')


%% Plot chi and eps profiles, for 5-profiles at a time, with all
% points and binned average profile

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

plot_dir = fullfile(fig_dir,'chi_eps_profiles_wPoints')
ChkMkDir(plot_dir)


Pmin=20;
screen_chi=1
screen_ml=1

for cnum=500:5:3000
    
    clear cnum_range
    cnum_range = [cnum cnum+5];
    
    clear cnums_to_get
    cnums_to_get = [cnum_range(1) : cnum_range(2) ];
    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);
    
    clear chipod cham
    [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
        path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
    %
    
    ib=find(log10(chipod.eps)>-4);
    chipod.eps(ib)=nan;
    ib=find(log10(cham.eps)>-4);
    cham.eps(ib)=nan;
    
    [chi_cham_avg, p_cham, ~] = binprofile(cham.chi, cham.P, 0, dz, 200,1);
    [chi_chi_avg, p_chi, ~] = binprofile(chipod.chi, chipod.P, 0, dz, 200,1);
    
    [eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
    [eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);
    
    
    figure(1);clf
    agutwocolumn(1)
    wysiwyg
    
    subplot(121)
    %h=scatter(log10(chipod.eps),chipod.P,'.','Markersize',15,'color','b')
    scatter(log10(chipod.chi),chipod.P,20,'b','filled','MarkerFaceAlpha',0.25)
    hold on
    %plot(log10(cham.eps),cham.P,'.','color',0.0*[1 1 1],'Markersize',15)
    h=scatter(log10(cham.chi),cham.P,20,'k','filled','MarkerFaceAlpha',0.25)
    hcham=plot(log10(chi_cham_avg),p_cham,'kd-','linewidth',2)
    hold on
    hchi=plot(log10(chi_chi_avg),p_chi,'bo-','linewidth',2)
    axis ij
    grid on
    legend([hcham hchi],'cham','\chi pod','location','best')
    %axis tight
    ylim([0 200])
    xlim([-11 -4])
    title(['profiles ' num2str(cnum_range(1)) ' - ' num2str(cnum_range(2))])
    xlabel('log_{10}[\chi]','fontsize',16)
    ylabel('P','fontsize',16)
    
    subplot(122)
    %h=scatter(log10(chipod.eps),chipod.P,'.','Markersize',15,'color','b')
    scatter(log10(chipod.eps),chipod.P,20,'b','filled','MarkerFaceAlpha',0.4)
    hold on
    %plot(log10(cham.eps),cham.P,'.','color',0.0*[1 1 1],'Markersize',15)
    h=scatter(log10(cham.eps),cham.P,20,'k','filled','MarkerFaceAlpha',0.4)
    hcham=plot(log10(eps_cham_avg),p_cham,'kd-','linewidth',2)
    hold on
    hchi=plot(log10(eps_chi_avg),p_chi,'bo-','linewidth',2)
    axis ij
    grid on
    legend([hcham hchi],'cham','\chi pod','location','best')
    %axis tight
    ylim([0 200])
    xlim([-8.5 -4])
    title(['profiles ' num2str(cnum_range(1)) ' - ' num2str(cnum_range(2))])
    xlabel('log_{10}[\epsilon]','fontsize',16)
    ylabel('P','fontsize',16)    
    
    print(fullfile(plot_dir,['eps_prof_cnum_' num2str(cnum)]),'-dpng')
    
    %pause(0.1)
    
end

%% calculate bias for 10m bins for 5-profile chunks, and try to figure
% out what is causing it (relation to N2,Tz,depth?)

clear ; close all

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=5; % bin size

eq14_patches_paths

plot_dir = fullfile(fig_dir,'chi_eps_profiles_wPoints')
ChkMkDir(plot_dir)

Pmin=20;
screen_chi=1
screen_ml=0

eps_cham = [];
eps_chi = [];
eps_bias = [];
P_all = [];
N2_all_chi = [];
N2_all_cham = [];
Tz_all = [];

for cnum=1:5:3000
    
    clear cnum_range
    cnum_range = [cnum cnum+1];
    
    clear cnums_to_get
    cnums_to_get = [cnum_range(1) : cnum_range(2) ];
    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);
    
    clear chipod cham
    [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
        path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
    %
    
    ib=find(log10(chipod.eps)>-5);
%    chipod.eps(ib)=nan;
    ib=find(log10(cham.eps)>-5);
%    cham.eps(ib)=nan;
    
    [chi_cham_avg, p_cham, ~] = binprofile(cham.chi, cham.P, 0, dz, 200,1);
    [chi_chi_avg, p_chi, ~] = binprofile(chipod.chi, chipod.P, 0, dz, 200,1);
    
    [eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
    [eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);

    [N2_cham_avg, p_cham, ~] = binprofile(cham.N2, cham.P, 0, dz, 200,1);
    [N2_chi_avg, p_chi, ~] = binprofile(chipod.N2, chipod.P, 0, dz, 200,1);

    [Tz_cham_avg, p_cham, ~] = binprofile(cham.Tz, cham.P, 0, dz, 200,1);
    [Tz_chi_avg, p_chi, ~] = binprofile(chipod.Tz, chipod.P, 0, dz, 200,1);
    
eps_chi = [eps_chi(:) ; eps_chi_avg(:)];
eps_cham = [eps_cham(:) ; eps_cham_avg(:)];
eps_bias = [ eps_bias(:) ; eps_chi_avg(:)./eps_cham_avg(:) ] ;
P_all = [ P_all(:) ; p_cham(:) ] ;
N2_all_chi = [N2_all_chi(:) ; N2_chi_avg(:)];
N2_all_cham = [N2_all_cham(:) ; N2_cham_avg(:)];
Tz_all = [Tz_all(:) ; Tz_chi_avg(:)];

end

%%
%id=find(log10(N2_all_chi)>-5);
id=find(P_all>80);
%
figure(1);clf
histogram(log10(eps_bias),'Normalization','pdf')
hold on
histogram(log10(eps_bias(id)),'Normalization','pdf')
grid on
xlim([-2 2])

%%

figure(22);clf
%loglog(N2_all_cham,N2_all_chi,'.')
%histogram2(log10(N2_all_cham), log10(N2_all_chi),'DisplayStyle','tile')
%histogram2(real(log10(N2_all_cham./N2_all_chi)),log10(eps_bias),'DisplayStyle','tile')
histogram2(real(log10(N2_all_cham./N2_all_chi)),P_all,'DisplayStyle','tile')
axis ij
%%

figure(2);clf
histogram2( log10(eps_cham), log10(eps_chi),25,'DisplayStyle','tile')
xlim([-9 -4])
ylim([-9 -4])
xvec=linspace(-9,-4,100);
hold on
plot(xvec,xvec,'k--')
xlabel('log_{10}\epsilon','fontsize',16)
ylabel('log_{10}\epsilon_{\chi}','fontsize',16)

%%

figure(3);clf
%plot(log10(eps_bias),P_all,'.')
histogram2(log10(eps_bias),P_all,'DisplayStyle','tile')
axis ij
xlim(2*[-1 1])
xlabel('eps bias')
ylabel('P')

%%

figure(4);clf
%plot(log10(eps_bias),P_all,'.')
histogram2(log10(eps_cham),log10(eps_bias),'DisplayStyle','tile')
%axis ij
ylim(2*[-1 1])
ylabel('eps bias')
xlabel('eps cham')

%%
figure(5);clf
%histogram2(log10(Tz_all),log10(eps_bias),'DisplayStyle','tile')
%histogram2(log10(N2_all_chi),log10(eps_bias),30,'DisplayStyle','tile')
histogram2(log10(N2_all_chi),P_all,'DisplayStyle','tile')
axis ij
%ylim(2*[-1 1])
ylabel('eps bias')
xlabel('N2')



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