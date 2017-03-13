
%%

clear ; close all

% load the 1m binned Cham data processed by Sally
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

iz=isin(cham.P,[60 200]);

figure(1);clf
scatter(cham.T1(:),cham.SAL(:),'.','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
hold on
scatter(cham.T1(iz),cham.SAL(iz),'.','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
ylim([34.2 35.6])
grid on

%%