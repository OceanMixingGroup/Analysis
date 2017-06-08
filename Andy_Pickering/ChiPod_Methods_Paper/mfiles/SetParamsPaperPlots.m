%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% SetParamsPaperPlots.m
%
% Set Params for all chi-pod paper plotting scripts. so we dont have to
% change parameters in every code...
%
%--------------------
% 07/10/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear Params cham

% load Chameleon data
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum.mat')load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum.mat')

Params.z_smooth=10;
Params.gamma=0.2
Params.fmax=7;
Params.resp_corr=0;
Params.fc=99;

%%