%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Processed_EQ14_cham.m
%
% Plot pre-processed (Sally/Jim) chameleon data from EQ14
%
%----------------------
% 01/19/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=1

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

m=4,n=1 ;

xl=[0 3100]
yl=[0 230]

ax1=subplot(m,n,1);
ezpc(cham.castnumber,cham.P,real(log10(cham.N2)))
colorbar
caxis([-5 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('N^2');

ax2=subplot(m,n,2);
ezpc(cham.castnumber,cham.P,real(log10(cham.DTDZ)))
colorbar
caxis([-3 -0])
xlim(xl)
ylim(yl)
SubplotLetterMW('dTdz');

ax3=subplot(m,n,3);
ezpc(cham.castnumber,cham.P,real(log10(cham.CHI)))
colorbar
caxis([-11 -3])
xlim(xl)
ylim(yl)
SubplotLetterMW('\chi');

ax4=subplot(m,n,4);
ezpc(cham.castnumber,cham.P,real(log10(cham.EPSILON)))
colorbar
caxis([-11 -4])
xlim(xl)
ylim(yl)
SubplotLetterMW('\epsilon');

linkaxes([ax1 ax2 ax3 ax4 ])

if saveplot==1
    %   figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures'
    eq14_patches_paths
    print( fullfile( fig_dir,['EQ14_PreProc_Summary'] ) , '-dpng' )
end

%%