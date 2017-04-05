
%%


clear ; close all

%whN2dTdz = 'line'
%whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
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
for cnum=1:25:3000
    try
        h = PlotEpsProfileCompare_eq08(cnum,Params,patch_size_min,...
            usetemp,minR2,dz)
        
        print( fullfile( figdir2, ['eq14_profile_' num2str(cnum) '_eps_profiiles_compare'] ),'-dpng')
        %pause(1)
    catch
        disp(['error on profile ' num2str(cnum)])
    end
    
end


%% compute 10m avg profiles and plot ratio of chameleon to binned profiles

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=32

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz = 10 % bin size

eq08_patches_paths

path_chipod_patches = fullfile(analysis_dir,project_long,'data','ChipodPatches');
path_chipod_bin     = '/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/';
path_cham_avg       = '/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/avg';

binned = [] ;
cham   = [] ;
patch  = [] ;

for cnum = 1:3000
    
    clear avg ch chb
    
    try
        % patch N^2,dTdz w/ constant gamma
        %         clear avg
        %         load( fullfile( path_chipod_patches, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        %         ch = avg;clear avg
        
        % regular chi-pod method on binned data
        clear avg
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        
        load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']))
        chb = avg;clear avg
        
        % chamelon data
        load(fullfile( path_cham_avg, ['eq08_' sprintf('%04d',cnum) '_avg.mat']) )
        
        clear bin1 bin2 bin3
        [bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
        %[bin3 z3 Nobs] = binprofile(ch.eps1    , ch.P , 0, dz, 200,0);
        
        cham   = [cham(:)   ; bin1(:) ];
        binned = [binned(:) ; bin2(:) ];
        %patch  = [ patch(:) ; bin3(:) ];
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum


% Histograms of chi-pod epsilon to chameleon epsilon (10bins)

figure(1);clf
h1 = histogram(log10(binned./cham),'EdgeColor','none','Normalization','pdf');
hold on
%h2 = histogram(log10(patch./cham),'EdgeColor','none','Normalization','pdf');
xlim([-4 4])
grid on
legend([h1 h2],'bin/cham','patch/cham')
xlabel('log_{10}[\epsilon ratio]')
ylabel('pdf')
title([project_short ' ' num2str(dz) ' m binned'])

%% plot average of groups of profiles, binned

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=32

% patch parameters
patch_size_min = 0.4
usetemp = 1
minR2 = 0.0

dz=10; % bin size

eq08_patches_paths


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
    cnums = [cnum_range(1) : cnum_range(2) ];
    
    %e1=[]; P1=[];
    ebin=[]; Pbin=[];
    echam=[]; Pcham=[];
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            %             clear avg
            %             load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
            %             e1 = [e1(:) ; avg.eps1(:)];
            %             P1 = [P1(:) ; avg.P(:)   ];
            
            % load binned chipod proifles
            clear avg
            path_chipod_bin = '/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/';
            load( fullfile( path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']))
            
            ebin = [ebin(:) ; avg.eps1(:)];
            Pbin = [Pbin(:) ; avg.P(:)   ];
            
            % load Chameleon profile
            clear avg
            path_cham_avg = '/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/avg/';
            load( fullfile(path_cham_avg,['eq08_' sprintf('%04d',cnum) '_avg.mat']))
            echam = [echam(:) ; avg.EPSILON(:)];
            Pcham = [Pcham(:) ; avg.P(:)   ];
            
        catch
        end
    end % cnum
    
    
    %    iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));
    
    ebin(find(log10(ebin)>-4))=nan;
    
    clear dataout1 dataout2 cham_bin
    %[dataout1 zout1 Nobs] = binprofile(e1,P1, 0, dz, 200,1);
    [ebin_avg zbin_avg Nobs] = binprofile(ebin,Pbin, 0, dz, 200,1);
    [cham_bin zout_cham Nobs] = binprofile(echam,Pcham, 0, dz, 200,1);
    
    %
    figure(1)
    subplot(2,3,whcase)
    %h1 = plot(log10(dataout1),zout1,'bo-','linewidth',2)
    %hold on
    h2 = plot(log10(ebin_avg),zbin_avg,'ms-','linewidth',2)
    hold on
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim([-9.5 -2])
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    legend([h2 h3],'bin','cham','location','best')
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    %eq14_patches_paths
    %print( fullfile(fig_dir,['eps_prof_cnums_' num2str(cnum_range(1)) '_' num2str(cnum_range(2))]),'-dpng')
    
    %     figure(2)
    %     subplot(3,2,whcase)
    %     loglog(cham_bin,dataout1,'o')
    %     hold on
    %     loglog(cham_bin,dataout2,'d')
    %     grid on
    %     xlim([1e-9 1e-5])
    %     ylim([1e-9 1e-5])
    %     xvec=linspace(1e-11,1e-4,100);
    %     plot(xvec,xvec,'k--')
    
    
end

%%
%figure(1)
%eq14_patches_paths
%print( fullfile(fig_dir,['eps_prof_comparisons']),'-dpng')


%%