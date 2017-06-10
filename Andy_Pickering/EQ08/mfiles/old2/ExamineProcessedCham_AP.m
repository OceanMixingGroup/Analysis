%%

clear ; close all

load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
%load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum_filled.mat')
%load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum_filtered.mat')

figure(1);clf
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
colorbar
caxis([-11 -4])

%%

%load('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/chi_all_zsm10m_fmax20Hz_respcorr1_gamma20.mat')
load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')

figure(1);clf

ax1=subplot(211)
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
colorbar
caxis([-11 -4])
ylim([0 210])

ax2=subplot(212)
ezpc(C.castnumber,C.p,log10(C.eps))
colorbar
caxis([-11 -4])
ylim([0 210])

linkaxes([ax1 ax2])

%% Plot a summary figure of the pre-processed Chameleon data

clear ; close all

saveplot=1

load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

m=4,n=1 ;

xl=[0 2680]
yl=[0 210]

ax1=subplot(m,n,1);
ezpc(cham.castnumber,cham.P,real(log10(cham.N2)))
colorbar
caxis([-5 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('N^2');
title('EQ08 Chameleon Summary')
ylabel('P [db]','fontsize',16)

ax2=subplot(m,n,2);
ezpc(cham.castnumber,cham.P,real(log10(cham.DTDZ)))
colorbar
caxis([-3 -0])
xlim(xl)
ylim(yl)
SubplotLetterMW('dTdz');
ylabel('P [db]','fontsize',16)

ax3=subplot(m,n,3);
ezpc(cham.castnumber,cham.P,real(log10(cham.CHI)))
colorbar
caxis([-11 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('\chi');
ylabel('P [db]','fontsize',16)

ax4=subplot(m,n,4);
ezpc(cham.castnumber,cham.P,real(log10(cham.EPSILON)))
colorbar
caxis([-11 -4])
xlim(xl)
ylim(yl)
SubplotLetterMW('\epsilon');
xlabel('Castnumber','fontsize',16)
ylabel('P [db]','fontsize',16)

% ax5=subplot(m,n,5);
% ezpc(cham.castnumber,cham.P,real(log10(cham.KT)))
% colorbar
% caxis([-9 -2])
% ylim(yl)
% SubplotLetterMW('KT')

linkaxes([ax1 ax2 ax3 ax4 ])

if saveplot==1
   figdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures'
   print( fullfile( figdir,['EQ08_PreProc_Summary'] ) , '-dpng' )
end

%%