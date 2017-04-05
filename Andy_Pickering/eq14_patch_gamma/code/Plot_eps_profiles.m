%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_eps_profiles.m
%
% Plot profiles of epsilon from chameleon and chi-pod method. Compare
% time-average profiles for groups of casts?
%
%
%----------------
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot single profiles: chameleon, chi-pod binned and patch

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz=10

eq14_patches_paths
figdir2 = fullfile( fig_dir, 'eps_profiles', whN2dTdz);
ChkMkDir(figdir2)

for cnum=1:25:3000
    try
        h = PlotEpsProfileCompare_eq14(cnum,whN2dTdz,Params,patch_size_min,...
            usetemp,minR2,dz)
        %        print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_' whN2dTdz '_eps_profiiles_compare'] ),'-dpng')
        
        pause(1)
    end
end

%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz = 10 % bin size

eq14_patches_paths

path_chipod_patches = fullfile(analysis_dir,project_long,'data','ChipodPatches');
%path_chipod_bin     = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/';
%path_cham_avg       = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/';

binned = [] ;
cham   = [] ;
patch  = [] ;

for cnum = 1:3000
    
    clear avg ch chb
    
    try
        % patch N^2,dTdz w/ constant gamma
        clear avg
        load( fullfile( path_chipod_patches, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        ch = avg;clear avg
        
        % regular chi-pod method on binned data
        clear avg
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        
        load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        chb = avg;clear avg
        
        % chamelon data
        load(fullfile( path_cham_avg, ['eq14_' sprintf('%04d',cnum) '.mat']) )
        
        clear bin1 bin2 bin3
        [bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
        [bin3 z3 Nobs] = binprofile(ch.eps1    , ch.P , 0, dz, 200,0);
        
        cham   = [cham(:)   ; bin1(:) ];
        binned = [binned(:) ; bin2(:) ];
        patch  = [ patch(:) ; bin3(:) ];
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum


%% Histograms of chi-pod epsilon to chameleon epsilon (10bins)

figure(1);clf
h1 = histogram(log10(binned./cham),'EdgeColor','none','Normalization','pdf');
hold on
h2 = histogram(log10(patch./cham),'EdgeColor','none','Normalization','pdf');
xlim([-4 4])
grid on
legend([h1 h2],'bin/cham','patch/cham')
xlabel('log_{10}[\epsilon ratio]')
ylabel('pdf')
title(['eq14 ' num2str(dz) ' m binned'])

%%
figname=['eq14_' num2str(dz) 'mbinned_eps_ratios']
print( fullfile(fig_dir,figname),'-dpng')
SetNotesFigDir
print( fullfile(NotesFigDir,figname),'-dpng')




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

%% plot average of groups of profiles, binned

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz=10; % bin size

eq14_patches_paths

path_chipod_patches = fullfile(analysis_dir,project_long,'data','ChipodPatches');
%path_chipod_bin     = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/';
%path_cham_avg       = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/';

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])

figure(1);clf
agutwocolumn(1)
wysiwyg

figure(2);clf
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
    
    e1=[]; P1=[];
    e2=[]; P2=[];
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            
            % patch chipod profile
            clear avg
            %load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            load( fullfile( path_chipod_patches, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            
            e1 = [e1(:) ; avg.eps1(:)];
            P1 = [P1(:) ; avg.P(:)   ];
            
            % binned chipod profile
            clear avg
            %path_chipod_bin = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/'
            load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            e2 = [e2(:) ; avg.eps1(:)];
            P2 = [P2(:) ; avg.P(:)   ];
            
        end
    end
    
    
    iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));
    
    e2(find(log10(e2)>-4))=nan;
    
    clear dataout1 dataout2 cham_bin
    [dataout1 zout1 Nobs] = binprofile(e1,P1, 0, dz, 200,1);
    [dataout2 zout2 Nobs] = binprofile(e2,P2, 0, dz, 200,1);
    [cham_bin zout_cham Nobs] = binprofile(cham.EPSILON(:,iCham),cham.P(:,iCham), 0, dz, 200,1);
    
    %
    figure(1)
    subplot(2,3,whcase)
    h1 = plot(log10(dataout1),zout1,'bo-','linewidth',2)
    hold on
    h2 = plot(log10(dataout2),zout2,'ms-','linewidth',2)
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim([-9.5 -2])
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    legend([h1 h2 h3],'patch','bin','cham','location','best')
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    %eq14_patches_paths
    %print( fullfile(fig_dir,['eps_prof_cnums_' num2str(cnum_range(1)) '_' num2str(cnum_range(2))]),'-dpng')
    
    figure(2)
    subplot(3,2,whcase)
    loglog(cham_bin,dataout1,'o')
    hold on
    loglog(cham_bin,dataout2,'d')
    grid on
    xlim([1e-9 1e-5])
    ylim([1e-9 1e-5])
    xvec=linspace(1e-11,1e-4,100);
    plot(xvec,xvec,'k--')
    
    
end

%
figure(1)
eq14_patches_paths
print( fullfile(fig_dir,['eps_prof_comparisons']),'-dpng')



%% different version where I plot all points

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

eq14_patches_paths

%dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches');
path_chipod_patches = fullfile(analysis_dir,project_long,'data','ChipodPatches');
%path_chipod_bin     = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/';
%path_cham_avg       = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/';

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

for whcase=6%1:6
    
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
            cnum_range = [2500 2600];
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
            %load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            load( fullfile( path_chipod_patches, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            
            e1 = [e1(:) ; avg.eps1(:)];
            P1 = [P1(:) ; avg.P(:)   ];
            
            % binned chi-pod
            clear avg
            %            load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
            load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            
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
    
    %
    figure(1);clf
    %subplot(2,3,whcase)
    plot(log10(e2),P2,'.','color',0.75*[1 1 1])
    hold on
    h1 = plot(log10(dataout1),zout1,'bo-','linewidth',2)
    hold on
    h2 = plot(log10(dataout2),zout2,'ms-','linewidth',2)
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    legend([h1 h2 h3],'patch','bin','cham','location','best')
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    %eq14_patches_paths
    %print( fullfile(fig_dir,['eps_prof_cnums_' num2str(cnum_range(1)) '_' num2str(cnum_range(2))]),'-dpng')
    
    xl=[-11 1]
    figure(2);clf
    subplot(1,3,1)
    plot(log10(e2),P2,'.','color',0.75*[1 1 1])
    hold on
    h2 = plot(log10(dataout2),zout2,'ms-','linewidth',2)
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    legend([h1 h2 h3],'patch','bin','cham','location','best')
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    xlim(xl)
    
    subplot(1,3,2)
    plot(log10(e1),P1,'.','color',0.75*[1 1 1])
    hold on
    h1 = plot(log10(dataout1),zout1,'bo-','linewidth',2)
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim(xl)
    ylim([0 200])
    
    subplot(1,3,3)
    plot(log10(cham.EPSILON(:,iCham)),cham.P(:,iCham),'.','color',0.75*[1 1 1])
    hold on
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim(xl)
    ylim([0 200])
    
end

%%

print(fullfile(fig_dir,'profexample'),'-dpng')

%%



%% See if gamma computed from multi-profile averageds of N2,Tz,chi,eps is 0.2?
% compare to gamma computed from individual 1m data points in every profile

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz=10; % bin size

eq14_patches_paths

dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches');

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])

cnum_range = [2400 3000];

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
[Tz_bin  zout Nobs] = binprofile(cham.DTDZ_RHOORDER(:,iCham),cham.P(:,iCham), 0, dz, 200,1);
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


%%

figure(1);clf
plot(log10(cham_bin),zout_cham)
hold on
plot(log10(dataout1),zout1)
plot(log10(dataout2),zout2)
axis ij

%%

figure(2);clf
loglog(cham_bin,dataout2,'o')
hold on
grid on

%%
ig = find(P2>80 & log10(e2)>-4);
i2 = find(P2>80 & log10(e2)<-4);

%%

figure(1);clf
subplot(121)
histogram(log10(N2(i2)),'Normalization','pdf','EdgeColor','none')
hold on
histogram(log10(N2(ig)),'Normalization','pdf')
subplot(122)
histogram(log10(Tz(i2)),'Normalization','pdf','EdgeColor','none')
hold on
histogram(log10(Tz(ig)),'Normalization','pdf')

% looks like very large estimates of epsilon from binned chi-pod method are
% associated with larger N2 and smaller Tz

%%

figure(2);clf

ax1=subplot(121)
boxplot(log10(Tz(i2)))
ylabel('log_{10}[T_z]')
grid on
title('log_{10}[\epsilon] <-4')
%
ax2=subplot(122)
boxplot(log10(Tz(ig)))
grid on
ylabel('log_{10}[T_z]')
title('log_{10}[\epsilon] >-4')
%
linkaxes([ax1 ax2])
%%

figure(2);clf

ax1=subplot(121)
boxplot(log10(abs(fspd(:))))
ylabel('log_{10}[T_z]')
grid on


ax2=subplot(122)
boxplot(log10(abs(fspd(ig))))
grid on
ylabel('log_{10}[T_z]')
title('log_{10}[\epsilon] >-4')

linkaxes([ax1 ax2])


%% scatter N2 vs Tz, color by epsilon ?

figure(3);clf

ax1=subplot(121)
scatter(log10(N2(i2)),log10(Tz(i2)),'filled','MarkerFaceAlpha',0.1)
grid on
xlabel('N^2')
ylabel('Tz')
title('All data')

ax2=subplot(122)
scatter(log10(N2(ig)),log10(Tz(ig)),'filled','MarkerFaceAlpha',0.1)
grid on
xlabel('N^2')
ylabel('Tz')
title('log_{10}[\epsilon] >-2')

linkaxes([ax1 ax2])

%%

figure(32);clf
ax1=subplot(121)
scatter(log10(tpvar(i2)),log10(Tz(i2)),'filled','MarkerFaceAlpha',0.1)
ax2=subplot(122)
scatter(log10(tpvar(ig)),log10(Tz(ig)),'filled','MarkerFaceAlpha',0.1)
xlabel('TPvar')
ylabel('T_z')

linkaxes([ax1 ax2])

%%
figure(32);clf
ax1=subplot(121)
scatter(log10(tpvar(i2)),log10(N2(i2)),'filled','MarkerFaceAlpha',0.1)
ax2=subplot(122)
scatter(log10(tpvar(ig)),log10(N2(ig)),'filled','MarkerFaceAlpha',0.1)
xlabel('TPvar')
ylabel('N^2')

linkaxes([ax1 ax2])

%%

figure(3);clf

ax1=subplot(121)
boxplot(log10(tpvar(:)))
grid on

ax2=subplot(122)
boxplot(log10(tpvar(ig)))
grid on

linkaxes([ax1 ax2])



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
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])

figure(1);clf
agutwocolumn(1)
wysiwyg

e1_all = [] ;
e2_all = [] ;
cham_all = [] ;

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
    
    e1=[]; P1=[];
    e2=[]; P2=[]; N2=[] ; Tz=[] ; tpvar=[]; fspd=[]; cnall=[];
    
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            
            % patch chi-pod
            clear avg
            load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            e1 = [e1(:) ; avg.eps1(:)];
            P1 = [P1(:) ; avg.P(:)   ];
            
            % binned chi-pod
            clear avg
            load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
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
    
end

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