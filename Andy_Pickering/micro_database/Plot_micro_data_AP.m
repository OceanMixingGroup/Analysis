%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_micro_data_AP.m
%
% Modified from plotKtOverKr.m (by Amy Waterhouse)
%
%--------------------
% 5/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%
clear ; close all

my_data_dir = '/Users/Andy/Google Drive/ChiCalculations/data/'
fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/micro_database/figures'
%
addpath /Users/Andy/Cruises_Research/mixingsoftware/colormaps/
%
% Kt/krho (using eqs 1 +2)
% KT =1/2 * chi / < dT=dz >^2 (1)
% Kr =0.2 * eps / N^2 (2)

colz = cbrewer('qual','Paired',26);

for i = [1:3,5,10]
    
    close all
    clear tnm K presK Kt hrp hrp96 hrp97
    
    if i == 3
        
        clear eps chi N2 dTdz hrp
        
        tnm = 'Natre';
        %load('../../data/natre/natre_hrp_for_map.mat')
        load( fullfile( my_data_dir,'natre','natre_hrp_for_map.mat'))
        
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz = (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 1
        
        clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
        tnm = 'BBTRE (smooth)';
        %load('../../data/bbtre/bbtre_hrp_for_map.mat')
        load( fullfile( my_data_dir,'bbtre','bbtre_hrp_for_map.mat'))
        hrp96=hrp;  %[1:20,71:75]
        %load('../../data/bbtre/bbtre97_hrp_for_map.mat')
        load( fullfile( my_data_dir,'bbtre','bbtre97_hrp_for_map.mat'))
        %hrp97=hrp;  % 1:36
        % why isn't hrp97 used?
        
        % Krho?
        hrp_bbtre_smooth = 0.2*([hrp96.eps(:, [2:5,7:20,71:74])])./([hrp96.N2(:, [2:5,7:20,71:74])]);
        
        K = hrp_bbtre_smooth;
        presK=hrp96.pres(:, [2:5,7:20,71:74]);
        
        KT = 0.5 * hrp96.chit(1:end-1,[2:5,7:20,71:74]) ./ (diff(hrp96.temp(:,[2:5,7:20,71:74]))./diff(hrp96.pres(:,[2:5,7:20,71:74]))).^2;
        
        eps = hrp.eps(:, [2:5,7:20,71:74]);
        chi = hrp.chit(:, [2:5,7:20,71:74]);
        N2  = hrp.N2(:, [2:5,7:20,71:74]);
        dTdz= (diffs(hrp.temp(:, [2:5,7:20,71:74]))./diffs(hrp.pres(:, [2:5,7:20,71:74])));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
        
    elseif i == 2
        
        clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
        
        tnm = 'BBTRE (rough)';
        %load('../../data/bbtre/bbtre_hrp_for_map.mat')
        load( fullfile( my_data_dir,'bbtre','bbtre_hrp_for_map.mat'))
        hrp96=hrp;
        %load('../../data/bbtre/bbtre97_hrp_for_map.mat')
        load(fullfile(my_data_dir,'bbtre','bbtre97_hrp_for_map.mat'))
        hrp97=hrp;
        
        % Krho?
        % mixed hrp96 and hrp97?
        hrp_bbtre_rough = 0.2*([hrp96.eps(1:10815,[1,6,21:70,75]) hrp97.eps(:,1:end)])./([hrp96.N2(1:10815,[1,6,21:70,75]) hrp97.N2(:,1:end)]);
        K = hrp_bbtre_rough;
        
        presK = [hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)];
        
        KT = 0.5 * ([hrp96.chit(1:10815-1,[1,6,21:70,75]) hrp97.chit(1:end-1,1:end)])./  (diff([hrp96.temp(1:10815,[1,6,21:70,75]) hrp97.temp(:,1:end)])./diff([hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)])).^2;
        
        clear temp p
        temp = [hrp96.temp(1:10815,[1,6,21:70,75]) hrp97.temp(:,1:end)];
        p= [hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)];
        
        eps = [hrp96.eps(1:10815,[1,6,21:70,75]) hrp97.eps(:,1:end)];
        chi = [hrp96.chit(1:10815,[1,6,21:70,75]) hrp97.chit(:,1:end)];%hrp.chit;
        N2  = [hrp96.N2(1:10815,[1,6,21:70,75]) hrp97.N2(:,1:end)];%hrp.N2;
        dTdz= diffs(temp)./diffs(p);% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 4
        
        %** missing chi? **
        
        clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
        
        tnm = 'Ladder';
        %load('../../data/ladder/ladder_hrp_for_map.mat')
        load (fullfile( my_data_dir,'ladder','ladder_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
        
    elseif i == 5
        
        tnm = 'Graviluck';
        
        clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
        
        % only plot down to 2216m since below that one profile exists
        %load('../../data/graviluck/graviluck_hrp_for_map.mat')
        load(fullfile(my_data_dir,'graviluck','graviluck_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 6
        
        %** no chi **
        
        clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
        
        tnm= 'Fieberling';
        %load('../../data/fieberling/fieberling_hrp_for_map.mat')
        load( fullfile( my_data_dir,'fieberling','fieberling_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 7
        
        % ** no chi?
        
        clear hrp K KT eps chi N2 dTdz
        
        tnm= 'Dimes_{DP}';
        %load('../../data/dimes/dimes_hrp_for_map.mat')
        load( fullfile( my_data_dir,'dimes','dimes_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 8
        
        clear hrp K KT eps chi N2 dTdz
        
        tnm= 'Dimes_{West}';
        %load('../../data/dimes/dimes_pacific_hrp_for_map.mat')
        load( fullfile( my_data_dir,'dimes','dimes_pacific_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 9
        
        clear hrp K KT eps chi N2 dTdz
        
        tnm= 'Toto';
        %load('../../data/toto/toto_hrp_for_map.mat')
        load( fullfile( my_data_dir, 'toto','toto_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chit;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    elseif i == 10
        tnm = 'Geotraces';
        % * note chi, not chit
        clear K KT eps chi N2 dTdz eps_chi
        %load('../../data/geotraces/geotraces_hrp_for_map.mat');
        load( fullfile( my_data_dir,'geotraces','geotraces_hrp_for_map.mat'))
        K = 0.2 * hrp.eps ./ hrp.N2;
        presK = hrp.pres;
        
        KT = 0.5 *  hrp.chi(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
        
        eps = hrp.eps;
        chi = hrp.chi;
        N2  = hrp.N2;
        dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
        eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
        
    end
    
    %~~ make figures
    
    
    % Amy's original 2X2 figure
    %~~~~~~~~~~~
    figure('paperposition',[0 0 11 4],'color','w'); wysiwyg;
    ax=subplot(221);
    pcolor(1:size(presK,2),nanmean(presK,2)',log10(abs(K)));
    axis ij; colorbar;
    caxis([-8 -4])
    ylabel('Depth (m)');
    xlabel('log_{10} K_{\rho} (m^2 s^{-1})')
    grid on; shading flat;
    title(tnm)
    colormap(ax,cbrewer('seq','Blues',11))
    
    ax=subplot(222);
    pcolor(1:size(presK,2),nanmean(presK(1:end-1,:),2)',log10(KT));
    axis ij; colorbar;
    set(gca,'yticklabel',[])
    caxis([-8 -4])
    xlabel('log_{10} K_{T} (m^2 s^{-1})')
    grid on; shading flat;
    colormap(ax,cbrewer('seq','Blues',11))
    
    ax2=subplot(223);
    pcolor(1:size(presK,2),nanmean(presK(1:end-1,:),2)',log10(KT./abs(K(1:end-1,:))));
    axis ij; colorbar;
    caxis([-2 2])
    ylabel('Depth (m)')
    %set(gca,'yticklabel',[])
    xlabel('log_{10} K_{T} / K_{\rho} ')
    grid on; shading flat;
    colormap(ax2,parula)
    
    fullratio = log10(KT./abs(K(1:end-1,:)));
    fullratio = fullratio(:);
    edges = [-5:.2:5];
    [n]=histcounts(fullratio,edges);
    
    ax2=subplot(224);
    bar(edges(1:end-1)+diff(edges(1:2))/2,n);
    xlabel('log10(K_T/K_{\rho})')
    ylabel('number of counts');
    set(gca,'xtick',[-5:1:5])
    grid on;
    
    %export_fig('-dpng','-r100',['KToverK_' tnm '.png']);
    
    %~~~~~~~~~~~
    
    %%
    try
        %~~~~~~~~~~~
        figure(2);clf
        histogram(fullratio,'edgeColor','none','Normalization','pdf')
        grid on
        freqline(nanmean(fullratio),'k--')
        title(tnm)
        xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
        print( fullfile( fig_dir,['KToverK_hist_' tnm '.png']), '-dpng')
        xlim([-3 3])
        print( fullfile( fig_dir,['KToverK_hist_xlim' tnm '.png']), '-dpng')
        
    end
    %%
    
    try
        figure(3);clf
        histogram2( real(log10(eps)), real(log10(eps_chi)),'DisplayStyle','tile')
        xlim([-12 -5])
        ylim([-12 -5])
        xvec = linspace(-12,-5,100);
        hold on
        plot(xvec,xvec,'k--')
        plot(xvec,xvec-1,'r--')
        plot(xvec,xvec+1,'r--')
        title(tnm)
        xlabel('log_{10}[\epsilon]','fontsize',16)
        ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)
        print( fullfile( fig_dir,['eps_chi_VS_eps_' tnm '.png']), '-dpng')
        
    end
    %%
    addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
    
    try
        figure(4);clf
        h = chi_vs_eps_normalized_plot(eps, chi, N2, dTdz)
        title(tnm)
        print( fullfile( fig_dir,['chi_eps_Norm_' tnm '.png']), '-dpng')
        
    end
    
    
    %%
end
