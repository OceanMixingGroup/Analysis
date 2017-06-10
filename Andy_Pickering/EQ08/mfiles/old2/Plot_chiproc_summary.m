%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_chiproc_summary.m
%
%
%--------------------
% 03/16/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot a summary figure of the chi-pod method applied to chameleon data

clear ; close all

saveplot=1

Params.fmax=20;
Params.z_smooth=10;
Params.gamma=0.2;
Params.resp_corr=1;

%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chi_all_' num2str(z_smooth) 'm.mat'])
%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chi_all_' num2str(Params.z_smooth) 'm_fmax_' num2str(Params.fmax) '.mat'])

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/',...
   ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]));

%load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

m=4,n=1 ;

xl=[0 2680]
yl=[0 210]

ax1=subplot(m,n,1);
ezpc(C.castnumber,C.p,real(log10(C.N2)))
colorbar
caxis([-5 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('N^2');
title(['EQ08 chi-proc Summary - zsm=' num2str(Params.z_smooth) 'm, fmax=' num2str(Params.fmax) 'Hz' ])
ylabel('P [db]','fontsize',16)

ax2=subplot(m,n,2);
ezpc(C.castnumber,C.p,real(log10(C.dTdz)))
colorbar
caxis([-3 -0])
xlim(xl)
ylim(yl)
SubplotLetterMW('dTdz');
ylabel('P [db]','fontsize',16)

ax3=subplot(m,n,3);
ezpc(C.castnumber,C.p,real(log10(C.chi)))
colorbar
caxis([-11 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('\chi');
ylabel('P [db]','fontsize',16)

ax4=subplot(m,n,4);
ezpc(C.castnumber,C.p,real(log10(C.eps)))
colorbar
caxis([-11 -4])
xlim(xl)
ylim(yl)
SubplotLetterMW('\epsilon');
xlabel('Castnumber','fontsize',16)
ylabel('P [db]','fontsize',16)

% ax5=subplot(m,n,5);
% ezpc(C.castnumber,C.P,real(log10(C.KT)))
% colorbar
% caxis([-9 -2])
% ylim(yl)
% SubplotLetterMW('KT')

linkaxes([ax1 ax2 ax3 ax4 ])

if saveplot==1
   figdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures'
   print( fullfile( figdir,['EQ08_chiproc_Summary'] ) , '-dpng' )
end

%%
