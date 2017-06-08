%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_hist_chieps_eps_all.m
%
% Plot histogram of ratio of chi_epsilon / epsilon for different datasets
%
%----------------
% 5/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%


clear ; close all

eps_thresh = 1   % option to impose lower threshold on epsilon
eps_floor  = -10 % value to use ( discards log10(eps)<eps_floor )

% set directory paths once here
fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/micro_database/figures'

figure(2);clf
agutwocolumn(0.7)
wysiwyg

hs=[] ; % collect handles for legend

for i = [1:3,5,10]
    
    % close all
    clear tnm K presK Kt hrp hrp96 hrp97 eps_chi eps chi N2 dTdz
    
    clear mix
    mix = Load_micro_data_mat(i,eps_thresh,eps_floor) ;
    
    
    
    %~~ make figures
    
    %%
    %try
    %~~~~~~~~~~~
    figure(2)
    h=histogram( real( log10(mix.eps_chi ./ mix.eps) ),'Normalization','pdf','DisplayStyle','stair','Linewidth',2)
    hold on
    grid on
    xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
    xlim([-3 3])
    hold on
    %
    hs=[hs h] ;
    
end % i

%% add data from EQ14

%~~~

clear cham eps_chi

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

eps_chi = cham.N2 .* cham.CHI / 2 / 0.2 ./ (cham.DTDZ.^2) ;
if eps_thresh==1
    ib = find(log10(cham.EPSILON)<-8.5);
    %ib = find(cham.P<80);
    cham.EPSILON(ib)=nan;
    eps_chi(ib)=nan;
end

h=histogram(log10(eps_chi./cham.EPSILON),'Normalization','pdf','EdgeColor','none')

%~~~
hs=[hs h]


%~~~~ add IWISE11 data from Lou
clear cham hrp
clear N2 dTdz eps_chi TEM PRS LAT SAL CH1
load('/Users/Andy/Cruises_Research/ChiPod/IWISE11_And_vmp/Data/VMP/IWISE-ASS.mat')

% compute N^2

N2   = nan * ones(size(CH1)) ;
dTdz = nan * ones(size(CH1)) ;
for ip = 1:length(LAT)
    N2(:,ip) = [ sw_bfrq( SAL(:,ip),TEM(:,ip),PRS,LAT(ip) ) ; nan ] ;
    dTdz(:,ip) = diffs( TEM(:,ip)) ./ diffs(PRS) ;
end
%
eps_chi = N2 .* CH1 /2 /0.2 ./ (dTdz.^2);

if eps_thresh==1
    ib = find(log10(eps)<eps_floor);
    eps(ib)=nan;
    eps_chi(ib)=nan;
    chi(ib)=nan;
end


h=histogram( real(log10(eps_chi./EP1)),'Normalization','pdf','DisplayStyle','stair','Linewidth',2)%,'EdgeColor','none')

%~~~
hs=[hs h]


ylim([0 1])
ylabel('pdf')
legend(hs,'BBTRE (smooth)', 'BBTRE (rough)','Natre','Graviluck', 'Geotraces','EQ14','IWISE11')

if eps_thresh==1
    print( fullfile( fig_dir,['epschi_eps_hist_ALL_eps_thresh.png']), '-dpng')
else
    print( fullfile( fig_dir,['epschi_eps_hist_ALL.png']), '-dpng')
end

%%