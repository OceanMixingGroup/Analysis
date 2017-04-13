%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_TIWE_avg_comb.m
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

tiwe_patches_paths

load( fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/','tiwe_1mavg_combined.mat') )

addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(311);
ezpc(cham.cnum,cham.P,cham.T)
hold on
contour(cham.cnum,cham.P,cham.T,[14:2:28],'k')
xlabel('cnum')
ylabel('P')
colorbar
title('T')

ax2=subplot(312);
ezpc(cham.cnum,cham.P,cham.S)
hold on
contour(cham.cnum,cham.P,cham.S,[5],'k')
xlabel('cnum')
ylabel('P')
colorbar
title('S')

ax3=subplot(313);
ezpc(cham.cnum,cham.P,sw_pden(cham.S,cham.T,cham.P,0))
hold on
contour(cham.cnum,cham.P,sw_pden(cham.S,cham.T,cham.P,0),[1020:1026],'k')
xlabel('cnum')
ylabel('P')
colorbar
title('sgth')

linkaxes([ax1 ax2 ax3])

%fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
print(fullfile(fig_dir,'tiwe_avgCombine_T_S_sgth'),'-dpng')


%%

figure(2);clf
histogram2( log10(cham.CHI(:)), log10(cham.EPSILON(:)),'DisplayStyle','tile')
xlabel('log_{10}[\chi]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)
colorbar
caxis([0 1000])
title('tiwe 1m avg data')


%%

figure(3);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(411);
ezpc(cham.cnum,cham.P,real(log10(cham.N2)))
%xlabel('cnum')
ylabel('P')
colorbar
caxis([-6 -1])
title('log_{10}[N^2]')
%xtloff

ax2=subplot(412);
ezpc(cham.cnum,cham.P,real(log10(cham.DTDZ)))
%xlabel('cnum')
ylabel('P')
colorbar
caxis([-4 0])
title('log_{10}[dT/dz]')
%xtloff

ax3=subplot(413);
ezpc(cham.cnum,cham.P,log10(cham.CHI))
%xlabel('cnum')
ylabel('P')
colorbar
caxis([-11 -4])
title('log_{10}[\chi]')
%xtloff

ax4=subplot(414);
ezpc(cham.cnum,cham.P,log10(cham.EPSILON))
xlabel('cnum')
ylabel('P')
colorbar
caxis([-11 -5])
title('log_{10}[\epsilon]')

linkaxes([ax1 ax2 ax3 ax4])

tiwe_patches_paths
print(fullfile(fig_dir,'tiwe_avgCombine_N2_dtdz_chi_eps'),'-dpng')

SetNotesFigDir
print(fullfile(NotesFigDir,'tiwe_avgCombine_N2_dtdz_chi_eps'),'-dpng')

%%