%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_chi_eps_norm_all.m
%
%
%----------------
% 5/25/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all


eps_thresh=0
eps_floor=-10

my_data_dir = '/Users/Andy/Google Drive/ChiCalculations/data/'
fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/micro_database/figures'
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

figure(2);clf
agutwocolumn(1)
wysiwyg

hs=[] ; % collect handles for legend
iax=1
for i = [1:3,5,10]
    
    % close all
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
    
    
    if eps_thresh==1
        ib = find(log10(eps)<eps_floor);
        eps(ib)=nan;
        eps_chi(ib)=nan;
        chi(ib)=nan;
        % axlims=[eps_floor -6];
    else
        % axlims=[-12 -6];
        
    end
    
    %~~~~~~~~~~~
    
    
    
    figure(2)
    subplot(3,2,iax)
    h = chi_vs_eps_normalized_plot(eps, chi, N2, dTdz) ;
    title(tnm)
    
    if iax~=1
        legend('off')
    end
    
    iax = iax + 1;
    
end % i

% add EQ14 data

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')
%eps_chi = cham.N2 .* cham.CHI / 2 / 0.2 ./ (cham.DTDZ.^2) ;

if eps_thresh==1
    ib = find(log10(cham.EPSILON)<-8.5);
    %ib = find(cham.P<80);
    cham.EPSILON(ib)=nan;
%    eps_chi(ib)=nan;
end


subplot(3,2,6)
h = chi_vs_eps_normalized_plot(cham.EPSILON, cham.CHI, cham.N2, cham.DTDZ) ;
title('EQ14')
legend('off')
%%

if eps_thresh==1
    print( fullfile( fig_dir,['chi_eps_norm_ALL_eps_thresh.png']), '-dpng')
    
else
    print( fullfile( fig_dir,['chi_eps_norm_ALL.png']), '-dpng')
end
%%