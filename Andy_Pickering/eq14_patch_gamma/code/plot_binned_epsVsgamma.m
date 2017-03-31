%%


clear ; close all

% load the 1m binned Cham data processed by Sally
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
cham1 = cham; clear cham

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

gam1 = ComputeGamma(cham1.N2,cham1.DTDZ_RHOORDER,cham1.CHI,cham1.EPSILON);

ig=find(log10(cham1.EPSILON)>-8.5);
figure(1);clf
histogram2(log10(gam1(ig)),log10(cham1.EPSILON(ig)),'displaystyle','tile')
xlim([-3 0])
ylim([-8.5 -5])
hf=freqline(log10(0.2),'r--')
set(hf,'linewidth',2)
xlabel('log_{10}[\gamma]')
ylabel('log_{10}[\epsilon]')

%%

figure(1);clf
scatter(log10(gam1(ig)),log10(cham1.EPSILON(ig)),'filled','MarkerFaceAlpha',0.05)
%xlim([-3 0])
%ylim([-8.5 -5])
hf=freqline(log10(0.2),'r--')


%%