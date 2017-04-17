%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_eps_profiles_eq08.m
%
% This script contains several analyses/plots, mostly comparing epsilon
% profiles from chi-pod method to chameleon, for different averaging etc.
%
%
%------------
% 4/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%


clear ; close all

Params.gamma = 0.2;
Params.fmax=32

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

eq08_patches_paths
figdir2 = fullfile( fig_dir, 'eps_profiles');
ChkMkDir(figdir2)
dz=10

for cnum=1:50:3000
    try
        h = PlotEpsProfileCompare_eq08(cnum,Params,patch_size_min,...
            usetemp,minR2,dz)
        
        print( fullfile( figdir2, [project_short '_profile_' num2str(cnum) '_eps_profiiles_compare'] ),'-dpng')
        %pause(1)
    catch
        disp(['error on profile ' num2str(cnum)])
    end
    
end


%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32 ;

dz = 10 % bin size

eq08_patches_paths

cnums_to_get=1:3000;

[eps_cham, chi_cham, N2_cham, Tz_cham, eps_chi, chi_chi, N2_chi, Tz_chi] =...
    Get_binned_data(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short);
%%

% eps_chi  = [] ;
% eps_cham = [] ;
% 
% chi_chi = [];
% chi_cham = [];
% 
% N2_chi = [];
% N2_cham = [];
% 
% Tz_chi = [];
% Tz_cham = [];
% 
% 
% for cnum = 1:3000
%     
%     clear avg ch chb
%     
%     try
%         
%         % regular chi-pod method on binned data
%         clear avg
%         load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']))
%         chb = avg;clear avg
%         
%         % chamelon data (1m bins)
%         load(fullfile( path_cham_avg, ['eq08_' sprintf('%04d',cnum) '_avg.mat']) )
%         
%         clear bin1 bin2
%         [bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
%         [bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
%         
%         eps_cham = [eps_cham(:)   ; bin1(:) ];
%         eps_chi  = [eps_chi(:)    ; bin2(:) ];
%         
%         clear bin1 bin2
%         [bin1 z1 Nobs] = binprofile(avg.CHI, avg.P, 0, dz, 200,1);
%         [bin2 z2 Nobs] = binprofile(chb.chi1   , chb.P, 0, dz, 200,1);
%         
%         chi_cham = [chi_cham(:)   ; bin1(:) ];
%         chi_chi  = [chi_chi(:)    ; bin2(:) ];
%         
%         clear bin1 bin2
%         [bin1 z1 Nobs] = binprofile(avg.N2, avg.P, 0, dz, 200,1);
%         [bin2 z2 Nobs] = binprofile(chb.N2   , chb.P, 0, dz, 200,1);
%         
%         N2_cham = [N2_cham(:)   ; bin1(:) ];
%         N2_chi  = [N2_chi(:)    ; bin2(:) ];
%         
%         clear bin1 bin2
%         [bin1 z1 Nobs] = binprofile(avg.DTDZ, avg.P, 0, dz, 200,1);
%         [bin2 z2 Nobs] = binprofile(chb.dTdz   , chb.P, 0, dz, 200,1);
%         
%         Tz_cham = [Tz_cham(:)   ; bin1(:) ];
%         Tz_chi  = [Tz_chi(:)    ; bin2(:) ];
%         
%     catch
%         disp(['error on profile ' num2str(cnum) ])
%     end % try
%     
% end % cnum
% 

%% Plot histogram of ratio of chi-pod epsilon to chameleon epsilon (10bins)

figure(1);clf
h1 = histogram(log10(eps_chi./eps_cham),'EdgeColor','none','Normalization','pdf');
hold on
xlim([-4 4])
grid on
%legend([h1 h2],'bin/cham','patch/cham')
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
ylabel('pdf')
title([project_short ' ' num2str(dz) ' m binned'])

eq08_patches_paths
print( fullfile(fig_dir,[project_short '_' num2str(dz) 'mbinned_eps_ratios']),'-dpng')

%% Plot eps vs chi, normalized so slope is equal to 1/2*gamma

ib = find(log10(eps_cham)<-8.5);
eps_cham(ib)=nan;

% Nan values in mixed layer
Pmin = 80;
ib=find(P_cham_avg<Pmin);
eps_cham(ib)=nan;

figure(2);clf
agutwocolumn(0.8)
wysiwyg

hh=histogram2( real(log10(eps_cham./N2_cham)), log10(chi_cham./(Tz_cham.^2)), 80, 'DisplayStyle','tile')
%hh=histogram2( log10(chi_chi./(Tz_chi.^2)) , real(log10(eps_chi./N2_chi)),80,'DisplayStyle','tile')
%hh=scatter( log10(chi_cham./(Tz_cham.^2)) , log10(eps_cham./N2_cham),'filled','MarkerFaceAlpha',0.05)
%loglog(chi_cham./(Tz_cham.^2),eps_cham./N2_cham,'.','color',0.75*[1 1 1])
grid on
hold on
xvec=linspace(1e-7,1e0,100);
h1=plot( log10(xvec), log10(xvec*2*0.2),'k-')
h2=plot( log10(xvec), log10(xvec*2*0.1),'r-')
h3=plot( log10(xvec), log10(xvec*2*0.05),'c-')
ylim([-7.5 0.5])
xlim([-5.5 0.5])
ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
title([project_short ' 10m binned '])

print( fullfile(fig_dir,[project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized']),'-dpng')
%%
%%

%% Now want to do same, but average some profiles together to see if 
% converges to gamma=0.2

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 32  ;

dz = 10 % bin size

Pmin=80

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg
iax=1;
for dp=[10 50 100 500]
    
    normx_all = [];
    normy_all = [];
    
    for ix = 1:round(3000/dp)%30
        
        
        clear cnums_to_get
        cnums_to_get = [ (ix-1)*dp : (ix*dp) ] ;
        
        %clear eps_cham chi_cham N2_cham Tz_cham
        %clear eps_chi chi_chi N2_chi Tz_chi
        %    [eps_cham, chi_cham, N2_cham, Tz_cham, eps_chi, chi_chi, N2_chi, Tz_chi] =...
        %        Get_binned_data(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short);
        
        clear eps_cham_avg chi_cham_avg N2_cham_avg Tz_cham_avg
        clear eps_chi_avg chi_chi_avg N2_chi_avg Tz_chi_avg
        [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg] =...
            Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin);
                
        normy_all = [ normy_all(:) ; chi_cham_avg./(Tz_cham_avg.^2)] ;
        normx_all = [ normx_all(:) ; eps_cham_avg./N2_cham_avg ];
        
    end % idx
    
    subplot(2,2, iax)    
    %hh=histogram2( log10(normx_all) , real(log10(normy_all)),80,'DisplayStyle','tile')
    %hh=histogram2( log10(chi_chi./(Tz_chi.^2)) , real(log10(eps_chi./N2_chi)),80,'DisplayStyle','tile')
    hh=scatter( log10(normx_all) , log10(normy_all),'filled','MarkerFaceAlpha',0.1)
    %loglog(chi_cham./(Tz_cham.^2),eps_cham./N2_cham,'.','color',0.75*[1 1 1])
    grid on
    hold on
    xvec=linspace(1e-7,1e0,100);
    h1=plot( log10(xvec), log10(xvec*2*0.2) ,'k-');
    h2=plot( log10(xvec), log10(xvec*2*0.1) ,'r-');
    h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
    ylim([-7.5 -1])
    xlim([-5.5 -1])
    ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
    xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
    legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
    title([project_short ' 10m binned ,' num2str(dp) ' profile averages'])
    
    %     fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_' num2str(dp) 'profileAvg']
    %     print( fullfile(fig_dir,fname),'-dpng')
    
    iax=iax+1;
    
end % dp (# profiles averaged together)
%
fname = [project_short '_' num2str(dz) 'mbinned_eps_vs_chi_normalized_diffNprofileAvg']
print( fullfile(fig_dir,fname),'-dpng')



%%
%% plot average of groups of profiles, binned

% 1st we get all data and pressures for a set of profiles
% then we average all the data in 10m depth bins

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32 ;

dz=10; % bin size

eq08_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg

figure(2);clf
agutwocolumn(1)
wysiwyg

figure(3);clf
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
    cnums = [cnum_range(1) : cnum_range(2) ];
    
    eps_chi_all=[]; P_chi_all=[];
    eps_cham_all=[]; P_cham_all=[];
    
    chi_chi_all=[];
    chi_cham_all=[];
    
    N2_chi_all=[];
    N2_cham_all=[];
    
    Tz_chi_all=[];
    Tz_cham_all=[];
    
    clear Ngood
    Ngood=0;
    
    for i=1:length(cnums)
        
        try
            cnum=cnums(i);
            
            % load binned chipod proifles
            clear avg
            load( fullfile( path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']))
            
            eps_chi_all = [eps_chi_all(:) ; avg.eps1(:)];
            chi_chi_all = [chi_chi_all(:) ; avg.chi1(:)];
            N2_chi_all  = [N2_chi_all(:)  ; avg.N2(:)  ];
            Tz_chi_all  = [Tz_chi_all(:)  ; avg.dTdz(:)];
            P_chi_all   = [P_chi_all(:)   ; avg.P(:)   ];
            
            % load Chameleon profile
            clear avg
            load( fullfile(path_cham_avg,['eq08_' sprintf('%04d',cnum) '_avg.mat']))
            eps_cham_all = [eps_cham_all(:) ; avg.EPSILON(:)];
            chi_cham_all = [chi_cham_all(:) ; avg.CHI(:)];
            N2_cham_all  = [N2_cham_all(:)  ; avg.N2(:)];
            Tz_cham_all  = [Tz_cham_all(:)  ; avg.DTDZ(:)];
            P_cham_all = [P_cham_all(:) ; avg.P(:)   ];
            
            Ngood=Ngood+1;
        catch
        end
    end % cnum
    
    
    eps_cham_all(find(log10(eps_cham_all)<-8.5))=nan;
    
    eps_chi_all(find(log10(eps_chi_all)>-4))=nan;
    
    clear eps_chi pbin_chi cham_bin chi_chi N2_chi Tz_chi
    [eps_chi pbin_chi Nobs] = binprofile(eps_chi_all,P_chi_all, 0, dz, 200,1);
    [chi_chi pbin_chi Nobs] = binprofile(chi_chi_all,P_chi_all, 0, dz, 200,1);
    [N2_chi  pbin_chi Nobs] = binprofile(N2_chi_all ,P_chi_all, 0, dz, 200,1);
    [Tz_chi  pbin_chi Nobs] = binprofile(Tz_chi_all ,P_chi_all, 0, dz, 200,1);
    
    clear eps_cham chi_cham N2_cham Tz_cham
    [eps_cham pbin_cham Nobs] = binprofile(eps_cham_all,P_cham_all, 0, dz, 200,1);
    [chi_cham pbin_cham Nobs] = binprofile(chi_cham_all,P_cham_all, 0, dz, 200,1);
    [N2_cham  pbin_cham Nobs] = binprofile(N2_cham_all ,P_cham_all, 0, dz, 200,1);
    [Tz_cham  pbin_cham Nobs] = binprofile(Tz_cham_all ,P_cham_all, 0, dz, 200,1);
    %
    figure(1)
    
    subplot(2,3,whcase)
    h2 = plot(log10(eps_chi),pbin_chi,'ms-','linewidth',2)
    hold on
    h3 = plot(log10(eps_cham),pbin_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim([-9.5 -4])
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    if whcase==5
        legend([h2 h3],'bin','cham','location','best')
    end
    title([num2str(cnum_range(1)) '-' num2str(cnum_range(2)) ' ,Ngood=' num2str(Ngood)])
    
    %     figure(2)
    %     subplot(3,2,whcase)
    %     loglog(eps_cham,eps_chi,'o')
    %     hold on
    %     grid on
    %     xlim([1e-9 1e-5])
    %     ylim([1e-9 1e-5])
    %     xvec=linspace(1e-11,1e-4,100);
    %     plot(xvec,xvec,'k--')
    %     title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    %
    figure(3)
    subplot(3,2,whcase)
    loglog(chi_chi./(Tz_chi.^2),eps_chi./N2_chi,'o')
    hold on
    loglog(chi_cham./(Tz_cham.^2),eps_cham./N2_cham,'d')
    grid on
    xvec=linspace(1e-7,1e0,100);
    h1=plot(xvec,xvec/2/0.2,'k-');
    h2=plot(xvec,xvec/2/0.1,'r-');
    xlim([1e-7 1e0])
    ylim([1e-5 1e0])
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    
    
end

%%
figure(1)
eq08_patches_paths
print( fullfile(fig_dir,[project_short '_eps_prof_comparisons']),'-dpng')



%% Try varying # of profiles averaged to see how many it takes to converge

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 32  ;

eq08_patches_paths

%figure(1);clf
%agutwocolumn(1)
%wysiwyg

e2_all = [] ;
echam_all = [] ;

prof_start=700
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
    
    echam=[]; Pcham=[];
    e2=[]; P2=[]; %
    
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            
            % binned chi-pod
            clear avg
            load(fullfile(path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']))
            e2 = [e2(:) ; avg.eps1(:)];
            P2 = [P2(:) ; avg.P(:)   ];
            
            % chamelon data (1m bins)
            load(fullfile( path_cham_avg, ['eq08_' sprintf('%04d',cnum) '_avg.mat']) )
            echam = [echam(:) ; avg.EPSILON(:)];
            Pcham = [Pcham(:) ; avg.P(:)   ];
            
        catch
        end
    end % cnums
    
    echam(find(log10(echam)>-4))=nan;
    e2(find(log10(e2)>-4))=nan;
    
    clear dataout2 cham_bin
    [dataout2 zout2 Nobs] = binprofile(e2,P2, 0, 10, 200,1);
    [cham_bin zout_cham Nobs] = binprofile(echam,Pcham, 0, 10, 200,1);
    
    e2_all    = [e2_all dataout2] ;
    echam_all = [echam_all cham_bin ];
    
end % whcase

%

figure(1);clf
agutwocolumn(1)
wysiwyg

for iax=1:6
    
    subplot(2,3,iax)
    plot(log10(e2_all(:,iax)),zout2,'ro-')
    hold on
    plot(log10(echam_all(:,iax)),zout_cham,'k')
    axis ij
    grid on
    xlim([-11 -4])
title([num2str(prof_start) ' - ' num2str(prof_start+prof_extra(iax))])
end

print(fullfile(fig_dir,[project_short '_eps_prof_diffN']),'-dpng')

%%

%% See if gamma computed from multi-profile averageds of N2,Tz,chi,eps is 0.2?
% compare to gamma computed from individual 1m data points in every profile

clear ; close all

dz=10; % bin size

eq08_patches_paths

% load processed chameleon data
load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')

cnum_range = [1400 3000];

clear cnums
cnums = [cnum_range(1) : cnum_range(2) ];

eps = [];
chi = [];
N2  = [];
Tz  = [];

for i=1:length(cnums)
    
    %    try
    iCham=find( cham.castnumber==cnums(i) );
    if ~isempty(iCham)
        eps = [eps(:) ; cham.EPSILON(:,iCham) ] ;
        chi = [chi(:) ; cham.CHI(:,iCham) ] ;
        N2  = [N2(:) ; cham.N2(:,iCham) ] ;
        Tz  = [Tz(:) ; cham.DTDZ(:,iCham) ] ;
    end
    %    catch
    %    end
end

% compute gamma from these values
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
gam = ComputeGamma(N2,Tz,chi,eps);

% Now average profiles together in 10m bins and compute gamma from that

iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));

dz=10

[N2_bin  zout Nobs] = binprofile(cham.N2(:,iCham)           ,cham.P(:,iCham), 0, dz, 200,1);
[Tz_bin  zout Nobs] = binprofile(cham.DTDZ(:,iCham),cham.P(:,iCham), 0, dz, 200,1);
[chi_bin zout Nobs] = binprofile(cham.CHI(:,iCham)          ,cham.P(:,iCham), 0, dz, 200,1);
[eps_bin zout Nobs] = binprofile(cham.EPSILON(:,iCham)      ,cham.P(:,iCham), 0, dz, 200,1);

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
title('profile-averaged, 10mbin')

linkaxes([ax1 ax2])

print(fullfile(fig_dir,[project_short '_gamma_point_avg_box']),'-dpng')

%%


%%
