%
%
% Plot profiles and bootstrapped averages
%
%
%
%% Plot single profiles
% points and binned average profile

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

plot_dir = fullfile(fig_dir,'N2_Tz_chi_eps')
ChkMkDir(plot_dir)


Pmin=10;
screen_chi=1
screen_ml=1

for cnum = 1:25:3000
    
    try
        clear cnum_range
        %cnum_range = [cnum cnum+5];
        
        clear cnums_to_get
        cnums_to_get = cnum;%[cnum_range(1) : cnum_range(2) ];
        bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
        cnums_to_get = setdiff(cnums_to_get,bad_prof);
        
        clear chipod cham
        [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
            path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
        %
        
        ib=find(log10(chipod.eps)>-5);
        %chipod.eps(ib)=nan;
        ib=find(log10(cham.eps)>-5);
        %cham.eps(ib)=nan;
        
        [chi_cham_avg, p_cham, ~] = binprofile(cham.chi, cham.P, 0, dz, 200,1);
        [chi_chi_avg, p_chi, ~] = binprofile(chipod.chi, chipod.P, 0, dz, 200,1);
        
        [eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
        [eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);
        
        
        figure(1);clf
        agutwocolumn(1)
        wysiwyg
        
        subplot(221)
        hchi = plot(log10(chipod.N2),chipod.P)
        hold on
        hcham = plot(log10(cham.N2),cham.P)
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        xlabel('N^2')
        ylim([0 200])
        
        subplot(222)
        hchi = plot(log10(chipod.Tz),chipod.P)
        hold on
        hcham = plot(log10(cham.Tz),cham.P)
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        xlabel('Tz')
        ylim([0 200])
        
        subplot(223)
        hchi = plot(log10(chipod.eps),chipod.P,'b')
        hold on
        hcham = plot(log10(cham.eps),cham.P,'r')
        hchi = scatter(log10(chipod.eps),chipod.P,40,'b','filled','MarkerFaceAlpha',0.5)
        hcham = scatter(log10(cham.eps),cham.P,40,'filled','MarkerFaceAlpha',0.5)
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        xlabel('\epsilon')
        ylim([0 200])
        
        subplot(224)
        hcham = plot(log10(eps_cham_avg),p_cham,'r')
        hold on
        hchi = plot(log10(eps_chi_avg),p_chi,'b')
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        ylim([0 200])
        %
        print(fullfile(plot_dir,['eps_boot_prof_cnum_' num2str(cnum)]),'-dpng')
    catch
    end
    
end

%% Plot single bootstrapped profiles

% points and binned average profile

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

plot_dir = fullfile(fig_dir,'boot_eps_profiles')
ChkMkDir(plot_dir)


Pmin=10;
screen_chi=1
screen_ml=1

for cnum = 1:25:3000
    
    clear cnum_range
    
    clear cnums_to_get
    cnums_to_get = cnum ;
    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);
    
    clear chipod cham
    [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
        path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
    %
    
    ib=find(log10(chipod.eps)>-5);
    %chipod.eps(ib)=nan;
    ib=find(log10(cham.eps)>-5);
    %cham.eps(ib)=nan;
    
    figure(1);clf
    agutwocolumn(0.8)
    wysiwyg
    
    try
        %
        ax1 = subplot(121);
        hchi  = plot(log10(chipod.eps),chipod.P,'b');
        hold on
        hcham = plot(log10(cham.eps),cham.P,'r');
        hchi  = scatter(log10(chipod.eps),chipod.P,40,'b','filled','MarkerFaceAlpha',0.5);
        hcham = scatter(log10(cham.eps),cham.P,40,'filled','MarkerFaceAlpha',0.5);
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        xlabel('\epsilon')
        xlim([-9 -4])
        ylim([0 200])
        
        ax2 = subplot(122);
        [bb, zbins, ~] = BinAndBootProfile(cham.eps,cham.P,10,100,1);
        plot(log10(bb(:,2)),zbins,'rd','linewidth',4)
        hold on
        for iz = 1:length(zbins)
            line([log10(bb(iz,[1 3]))],[zbins(iz) zbins(iz)],'color','r')
        end
        
        [bb, zbins, ~] = BinAndBootProfile(chipod.eps,chipod.P,10,100,0);
        plot(log10(bb(:,2)),zbins+2,'bp','linewidth',4)
        hold on
        for iz = 1:length(zbins)
            line([log10(bb(iz,[1 3]))],[zbins(iz)+2 zbins(iz)+2],'color','b')
        end
        
        axis ij
        grid on
        ylabel('P')
        xlabel('\epsilon')
        title(['profile ' num2str(cnum)])
        xlim([-9 -4])
        ylim([0 200])
        
        print(fullfile(plot_dir,['eps_boot_prof_cnum_' num2str(cnum)]),'-dpng')
    catch
    end % try
    
end % cnum



%% Plot bootstrapped average of several profiles

% points and binned average profile

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

plot_dir = fullfile(fig_dir,'boot_eps_profiles_multprof')
ChkMkDir(plot_dir)


Pmin=10;
screen_chi=1
screen_ml=1

dp = 5
for cnum = 1:25:3000
    
    clear cnum_range
    cnum_range = [cnum cnum+dp];
    
    clear cnums_to_get
    cnums_to_get = [cnum_range(1) : cnum_range(2) ];
    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);
    
    clear chipod cham
    [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
        path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
    %
    
    ib=find(log10(chipod.eps)>-5);
    %chipod.eps(ib)=nan;
    ib=find(log10(cham.eps)>-5);
    %cham.eps(ib)=nan;
    
    
    figure(1);clf
    agutwocolumn(0.8)
    wysiwyg
    
    try
        %
        ax1 = subplot(121) ;
        hchi = scatter(log10(chipod.eps),chipod.P,40,'b','filled','MarkerFaceAlpha',0.5)
        hold on
        hcham = scatter(log10(cham.eps),cham.P,40,'filled','MarkerFaceAlpha',0.5)
        axis ij
        grid on
        legend([hchi hcham],'\chipod','cham','location','best')
        xlabel('\epsilon')
        xlim([-9 -4])
        ylim([0 200])
        
        ax2 = subplot(122);
        [bb, zbins, ~] = BinAndBootProfile(cham.eps,cham.P,10,100,1);
        plot(log10(bb(:,2)),zbins,'rd','linewidth',4)
        hold on
        for iz = 1:length(zbins)
            line([log10(bb(iz,[1 3]))],[zbins(iz) zbins(iz)],'color','r')
        end
        
        [bb, zbins, ~] = BinAndBootProfile(chipod.eps,chipod.P,10,100,0);
        plot(log10(bb(:,2)),zbins+2,'bp','linewidth',4)
        hold on
        %       plot(log10(bb(:,[1 3])),zbins,'rd','linewidth',4)
        for iz = 1:length(zbins)
            line([log10(bb(iz,[1 3]))],[zbins(iz)+2 zbins(iz)+2],'color','b')
        end
        
        axis ij
        grid on
        ylabel('P')
        xlabel('\epsilon')
        title(['profile ' num2str(cnum) ' - ' num2str(cnum+dp)])
        xlim([-9 -4])
        ylim([0 200])
        
        linkaxes([ax1 ax2])
        print(fullfile(plot_dir,['eps_boot_prof_cnum_' num2str(cnum)]),'-dpng')
        
    catch
    end % try
    
end % cnum
%%