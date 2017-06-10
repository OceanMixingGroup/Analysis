%
%
%
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam = ComputeGamma(cham.N2,cham.DTDZ,cham.CHI,cham.EPSILON);

%%
figure(1);clf
histogram2(log10(gam(:)),cham.P(:),'DisplayStyle','tile')
xlim([-5 5])
freqline(log10(0.2))
freqline(log10(0.1))
axis ij

%%
figure(2);clf
histogram2( log10(cham.EPSILON(:)), log10(gam(:)),'DisplayStyle','tile')

%%