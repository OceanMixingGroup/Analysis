%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakeBiasTable.m
%
% Make a table with the chipod/cham bias in different epsilon ranges
%
%-----------
% 08/05/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

% Set chipod params and load chameleon data
SetParamsPaperPlots
Params
%Params.nfft=256
Params.nfft=128

%~ load chi-pod processed data
%load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
%    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft)]));

% average chameleon data in same depth bins as chipod method
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)

%%

%ib=find(C.p<10);
%C.chi(ib,:)=nan;

chirat= C.chi ./ cham_bin.chi ;
LC=log10(chirat);
LE=log10(cham_bin.eps);

%%

epsranges=[-9 -8 ; -8 -7 ; -7 -6 ; -6 -5 ]
and=' & '
lend=' \\ '
clc

disp(['$\epsilon_{cham}$ range' and 'bias' and 'sd' lend])
disp('\hline')
for i=1:size(epsranges,1)
clear ig
ig=find( LE>epsranges(i,1) & LE<epsranges(i,2));
disp(['$' num2str(epsranges(i,1)) ' < \epsilon < ' num2str(epsranges(i,2)) '$' and num2str(roundx(nanmean(LC(ig)),2)) and num2str(roundx(nanstd(LC(ig)),2)) lend ])
disp('\hline')
end

%%