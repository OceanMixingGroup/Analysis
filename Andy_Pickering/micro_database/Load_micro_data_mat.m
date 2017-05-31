function mix = Load_micro_data_mat(i,eps_thresh,eps_floor)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% mix = Load_micro_data_mat(i,eps_thresh,eps_floor)
%
% Function to load data from Amy's microstructure database. This loads mat
% data that Amy gave me, not the raw .nc data from the database.
%
% From part of Amy's original script to make plots of KT/K etc..
%
% INPUT
% - i
% - eps_thresh
% - eps_floor
%
% OUTPUT
% - mix : Structure with data
%
%
%---------------
% 5/31/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%for i = [1:3,5,10]

my_data_dir = '/Users/Andy/Google Drive/ChiCalculations/data/'


% close all
clear tnm K presK Kt hrp hrp96 hrp97 eps_chi eps chi N2 dTdz


if i == 1
    
    clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
    tnm = 'BBTRE (smooth)';
    %load('../../data/bbtre/bbtre_hrp_for_map.mat')
    load( fullfile( my_data_dir,'bbtre','bbtre_hrp_for_map.mat'))
    hrp96 = hrp;  %[1:20,71:75]
    
    % Krho?
    %hrp_bbtre_smooth = 0.2 * hrp96.eps(:, [2:5,7:20,71:74]) ./ hrp96.N2(:, [2:5,7:20,71:74] );
    
    %K = hrp_bbtre_smooth;
    %presK=hrp96.pres(:, [2:5,7:20,71:74]);
    
    %KT = 0.5 * hrp96.chit(1:end-1,[2:5,7:20,71:74]) ./ (diff(hrp96.temp(:,[2:5,7:20,71:74]))./diff(hrp96.pres(:,[2:5,7:20,71:74]))).^2;
    
    eps = hrp.eps(:, [2:5,7:20,71:74]);
    chi = hrp.chit(:, [2:5,7:20,71:74]);
    N2  = hrp.N2(:, [2:5,7:20,71:74]);
    dTdz= (diffs(hrp.temp(:, [2:5,7:20,71:74]))./diffs(hrp.pres(:, [2:5,7:20,71:74])));
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
    
elseif i == 2
    
    clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
    
    tnm = 'BBTRE (rough)';
    %load('../../data/bbtre/bbtre_hrp_for_map.mat')
    load( fullfile( my_data_dir,'bbtre','bbtre_hrp_for_map.mat'))
    hrp96 = hrp;
    %load('../../data/bbtre/bbtre97_hrp_for_map.mat')
    load(fullfile(my_data_dir,'bbtre','bbtre97_hrp_for_map.mat'))
    hrp97 = hrp;
    
    clear temp p
    temp = [hrp96.temp(1:10815,[1,6,21:70,75]) hrp97.temp(:,1:end)];
    p= [hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)];
    
    eps = [ hrp96.eps(1:10815,[1,6,21:70,75]) hrp97.eps(:,1:end) ];
    chi = [ hrp96.chit(1:10815,[1,6,21:70,75]) hrp97.chit(:,1:end) ];
    N2  = [ hrp96.N2(1:10815,[1,6,21:70,75]) hrp97.N2(:,1:end) ];
    dTdz= diffs(temp)./diffs(p);
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
    % Krho?
    % mixed hrp96 and hrp97?
    %hrp_bbtre_rough = 0.2*([hrp96.eps(1:10815,[1,6,21:70,75]) hrp97.eps(:,1:end)])./([hrp96.N2(1:10815,[1,6,21:70,75]) hrp97.N2(:,1:end)]);
    %K = hrp_bbtre_rough;
    K = 0.2 * eps ./ N2 ;
    
    presK = [hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)];
    
    %KT = 0.5 * ([hrp96.chit(1:10815-1,[1,6,21:70,75]) hrp97.chit(1:end-1,1:end)])./  (diff([hrp96.temp(1:10815,[1,6,21:70,75]) hrp97.temp(:,1:end)])./diff([hrp96.pres(1:10815,[1,6,21:70,75]) hrp97.pres(:,1:end)])).^2;
    KT = 0.5 * chi ./ (dTdz.^2) ;
    
    
elseif i == 3
    
    clear eps chi N2 dTdz hrp
    
    tnm = 'Natre';
    %load('../../data/natre/natre_hrp_for_map.mat')
    load( fullfile( my_data_dir,'natre','natre_hrp_for_map.mat'))
    
    %         K = 0.2 * hrp.eps ./ hrp.N2;
    %         presK = hrp.pres;
    %
    %         KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    eps  = hrp.eps ;
    chi  = hrp.chit;
    N2   = hrp.N2  ;
    dTdz = (diffs(hrp.temp)./diffs(hrp.pres)) ;
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2) ;
    
    
elseif i == 4
    
    %** missing chi? **
    
    clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
    
    tnm = 'Ladder';
    %load('../../data/ladder/ladder_hrp_for_map.mat')
    load (fullfile( my_data_dir,'ladder','ladder_hrp_for_map.mat'))
    
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    eps = hrp.eps ;
    chi = hrp.chit;
    N2  = hrp.N2  ;
    dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
    
elseif i == 5
    
    tnm = 'Graviluck';
    
    clear hrp hrp96 hrp97 K KT eps chi N2 dTdz
    
    % only plot down to 2216m since below that one profile exists
    %load('../../data/graviluck/graviluck_hrp_for_map.mat')
    load(fullfile(my_data_dir,'graviluck','graviluck_hrp_for_map.mat'))
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %        KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
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
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
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
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
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
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    eps = hrp.eps ;
    chi = hrp.chit;
    N2  = hrp.N2 ;
    dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
elseif i == 9
    
    clear hrp K KT eps chi N2 dTdz
    
    tnm= 'Toto';
    %load('../../data/toto/toto_hrp_for_map.mat')
    load( fullfile( my_data_dir, 'toto','toto_hrp_for_map.mat'))
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chit(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    eps = hrp.eps ;
    chi = hrp.chit;
    N2  = hrp.N2  ;
    dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
elseif i == 10
    tnm = 'Geotraces';
    % * note chi, not chit
    clear K KT eps chi N2 dTdz eps_chi
    %load('../../data/geotraces/geotraces_hrp_for_map.mat');
    load( fullfile( my_data_dir,'geotraces','geotraces_hrp_for_map.mat'))
    
    %K = 0.2 * hrp.eps ./ hrp.N2;
    %presK = hrp.pres;
    
    %KT = 0.5 *  hrp.chi(1:end-1,:) ./ (diff(hrp.temp)./diff(hrp.pres)).^2;
    
    eps = hrp.eps;
    chi = hrp.chi;
    N2  = hrp.N2 ;
    dTdz= (diffs(hrp.temp)./diffs(hrp.pres));% dTdz = [dTdz(:) ; nan];
    eps_chi = N2 .* chi / 2 / 0.2 ./ (dTdz.^2);
    
end


%
if eps_thresh==1
    ib = find(log10(eps)<eps_floor);
    eps(ib)=nan;
    eps_chi(ib)=nan;
    chi(ib)=nan;
end

mix = struct('chi',chi,'eps',eps,'eps_chi',eps_chi,'N2',N2,'dTdz',dTdz,'project',tnm) ;

return

