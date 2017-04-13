%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_eps_profiles_tiwe.m
%
%
%
%------------
% 4/7/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32

dz=10 % dz for binning

tiwe_patches_paths

figdir2 = fullfile( fig_dir, 'eps_profiles');
ChkMkDir(figdir2)

for cnum=1:50:3000
    try
        h = PlotEpsProfileCompare_tiwe(cnum,Params,dz)
        
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

tiwe_patches_paths

path_chipod_patches = fullfile(analysis_dir,project_long,'data','ChipodPatches');

binned = [] ;
cham   = [] ;
patch  = [] ;

for cnum = 1:3000
    
    clear avg ch chb
    
    try
        
        % regular chi-pod method on binned data
        clear avg
        load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
        chb = avg;clear avg
        
        % chamelon data (1m bins)
        load(fullfile( path_cham_avg, ['tw91_' sprintf('%04d',cnum) '_avg.mat']) )
        
        clear bin1 bin2 
        [bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
        
        cham   = [cham(:)   ; bin1(:) ];
        binned = [binned(:) ; bin2(:) ];
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum


% Plot histogram of ratio of chi-pod epsilon to chameleon epsilon (10bins)

figure(1);clf
h1 = histogram(log10(binned./cham),'EdgeColor','none','Normalization','pdf');
hold on
xlim([-4 4])
grid on
%legend([h1 h2],'bin/cham','patch/cham')
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
ylabel('pdf')
title([project_short ' ' num2str(dz) ' m binned'])

%% plot average of groups of profiles, binned

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 32 ;

dz=10; % bin size

tiwe_patches_paths

figure(1);clf
agutwocolumn(1)
wysiwyg

figure(2);clf
agutwocolumn(1)
wysiwyg


for whcase=1:8
    
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
        case 7
            cnum_range = [3000 3500];
        case 8
            cnum_range = [3500 4000];
    end
    
    clear cnums
    cnums = [cnum_range(1) : cnum_range(2) ];
    
    ebin=[]; Pbin=[];
    echam=[]; Pcham=[];
    
    for i=1:length(cnums)
        try
            cnum=cnums(i);
            
            % load binned chipod proifles
            clear avg
            load( fullfile( path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
            
            ebin = [ebin(:) ; avg.eps1(:)];
            Pbin = [Pbin(:) ; avg.P(:)   ];
            
            % load Chameleon profile
            clear avg
            load( fullfile(path_cham_avg,['tw91_' sprintf('%04d',cnum) '_avg.mat']))
            echam = [echam(:) ; avg.EPSILON(:)];
            Pcham = [Pcham(:) ; avg.P(:)   ];
            
        catch
        end
    end % cnum
    
        
    ebin(find(log10(ebin)>-4))=nan;
    
    clear dataout1 dataout2 cham_bin
    [ebin_avg zbin_avg Nobs ] = binprofile(ebin,Pbin, 0, dz, 200,1);
    [cham_bin zout_cham Nobs] = binprofile(echam,Pcham, 0, dz, 200,1);
    
    %
    figure(1)
    
    subplot(4,2,whcase)
    h2 = plot(log10(ebin_avg),zbin_avg,'ms-','linewidth',2)
    hold on
    h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
    axis ij
    grid on
    xlim([-9.5 -4])
    ylim([0 200])
    xlabel('log_{10}[\epsilon]')
    ylabel('P [db]')
    if whcase==5
    legend([h2 h3],'bin','cham','location','best')
    end
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    
    figure(2)
    subplot(4,2,whcase)
    loglog(cham_bin,ebin_avg,'o')
    hold on
    grid on
    xlim([1e-9 1e-5])
    ylim([1e-9 1e-5])
    xvec=linspace(1e-11,1e-4,100);
    plot(xvec,xvec,'k--')
    title(['cnums ' num2str(cnum_range(1)) '-' num2str(cnum_range(2))])
    
    
end

%%

figure(1)
tiwe_patches_paths
print( fullfile(fig_dir,[project_short '_eps_prof_comparisons']),'-dpng')



%% Try varying # of profiles averaged to see how many it takes to converge

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 32  ;

tiwe_patches_paths

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
            load(fullfile(path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
            e2 = [e2(:) ; avg.eps1(:)];
            P2 = [P2(:) ; avg.P(:)   ];
            
            % chamelon data (1m bins)
            load(fullfile( path_cham_avg, ['tw91_' sprintf('%04d',cnum) '_avg.mat']) )
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
    title([num2str(prof_extra(iax))])
end

print(fullfile(fig_dir,[project_short 'eps_prof_diffN']),'-dpng')

%%