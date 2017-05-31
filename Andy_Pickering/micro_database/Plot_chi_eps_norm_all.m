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


eps_thresh =  1
eps_floor  = -10

%my_data_dir = '/Users/Andy/Google Drive/ChiCalculations/data/'
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
    
    clear mix
    mix = Load_micro_data_mat(i,eps_thresh,eps_floor) ;
    
    
    %~~ make figures
    
    figure(2)
    subplot(3,2,iax)
    h = chi_vs_eps_normalized_plot(mix.eps, mix.chi, mix.N2, mix.dTdz) ;
    title(mix.project)
    
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