%%



%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam = ComputeGamma(cham.N2,cham.DTDZ,cham.CHI,cham.EPSILON);

figure(1);clf
histogram2(log10(gam(:)),log10(cham.EPSILON(:)),'DisplayStyle','tile')
freqline(log10(0.2))
xlim([-4 2])

%%
