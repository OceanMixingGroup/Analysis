%~~~~~~~~~~~~~~~~~~~~~
%
% PlotChiForDiffu.m
%
% Compare results for different ranges of u to show there isn't a bias
%
% Need to go back to CombineCasts and select ranges of u at that step
%
%----------
% 05/26/16 - AP
%~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

saveplot=0

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures'

Params.fmax=7;
Params.z_smooth=10;
Params.gamma=0.2;
Params.resp_corr=0;
Params.fc=99

% load chi-pod processed data
%load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
%    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]));

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));

%%