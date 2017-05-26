%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_iwise11_Lou_plots.m
%
%
%
%----------------
% 5/26/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/micro_database/figures'

load('/Users/Andy/Cruises_Research/ChiPod/IWISE11_And_vmp/Data/VMP/IWISE-ASS.mat')

% compute N^2

N2 = nan * ones(size(CH1));
dTdz =nan * ones(size(CH1));
for ip = 1:length(LAT)
    N2(:,ip) = [sw_bfrq(SAL(:,ip),TEM(:,ip),PRS,LAT(ip)) ; nan];
    dTdz(:,ip) = diffs( TEM(:,ip)) ./ diffs(PRS) ;
end
%
eps_chi = N2 .* CH1 /2 /0.2 ./ (dTdz.^2);

%%

figure(1);clf
histogram( real(log10(eps_chi./EP1)),'EdgeColor','None','Normalization','pdf')
grid on
xlim([-3 3])
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
title('IWISE 11 - ASS')

print( fullfile( fig_dir,['KToverK_hist_iwise11ASS.png']), '-dpng')
        
%%

figure(3);clf
histogram2( real(log10(EP1(:))), real(log10(eps_chi(:))),'DisplayStyle','tile')
xlim([-12 -5])
ylim([-12 -5])
xvec = linspace(-12,-5,100);
hold on
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
title('IWISE 11 - ASS')
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)

print( fullfile( fig_dir,['eps_chi_VS_eps_iwise11ASS.png']), '-dpng')
        
%%

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/


figure(4);clf
h = chi_vs_eps_normalized_plot(EP1, CH1, N2, dTdz)

print( fullfile( fig_dir,['chi_eps_Norm_iwise11ASS.png']), '-dpng')
     
%%