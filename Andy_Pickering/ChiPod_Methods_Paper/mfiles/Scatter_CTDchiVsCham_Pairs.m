%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Scatter_CTDchiVsCham_Pairs.m
%
% * Makes plot for chipod methods paper *
%
% Uses data made in Make_Combined_Cham_for_CTD_pairs.m
%
% Scatter-plot \chi from chipod versus the mean of bracketing Chameleon
% profiles, for EQ14 CTD-chipod casts.
%
%--------------
% 05/12/16 - A.Pickering
% 06/03/16 - AP - Clean up and comment
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% same as above, but add histogram to right

clear ; close all

saveplot=1

load('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Ctd_Cham_Pairs.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)
%subplot(1,3,[1 2])
subplot(211)
h3=loglog(Pairs.chi_mn(:),Pairs.CTD.chi1(:),'.');
grid on
hold on
xvec=linspace(1e-12,1e-3,100);
loglog(xvec,xvec,'k--')
loglog(xvec,xvec/10,'r--')
loglog(xvec,xvec*10,'r--')

xlim([1e-12 1e-3])
ylim([1e-12 1e-3])
xlabel('log_{10} \chi_{\epsilon}^{cham} [K^2s^{-1}]','fontsize',16)
ylabel('log_{10} \chi_{\chi}^{ctd}] [K^2s^{-1}] ','fontsize',16)
title('CTD-\chi pod vs. Cham. mean')
axis square
set(gca,'Xtick',[1e-12 1e-10 1e-8 1e-6 1e-4])
set(gca,'Ytick',[1e-12 1e-10 1e-8 1e-6 1e-4])
pos=get(gca,'Position');
set(gca,'Position',pos.*[1.1 1 0.9 1])

crat=Pairs.CTD.chi1(:) ./ Pairs.chi_mn(:) ;
%subplot(1,3,3)
subplot(212)
histogram( log10( crat),'edgecolor','none','Normalization','pdf');
grid on
xlim(4*[-1 1])
grid on
xlabel('log_{10}[\chi_{\chi}^{ctd} / \chi_{\epsilon}^{cham}]','fontsize',16)
ylabel('counts','fontsize',16)
freqline(nanmean(log10(crat)),'b--');
title(['\mu= ' num2str(roundx(nanmean(log10(crat)),2)) ' , \sigma=' num2str(roundx(nanstd(log10(crat)),2))])

if saveplot==1
    SetPaperFigPath
    figname=['EQ14_ctdChipod_vs_chamMean_chi_scatter_WITHhist']
    print( fullfile( figdir , figname ), '-dpng' )
end

%%
% %% same as above, but just plot mean of chameleon pairs
% 
% clear ; close all
% 
% saveplot=1
% 
% load('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Ctd_Cham_Pairs.mat')
% 
% figure(1);clf
% agutwocolumn(0.7)
% wysiwyg
% h3=loglog(Pairs.chi_mn(:),Pairs.CTD.chi1(:),'.')
% grid on
% hold on
% xvec=linspace(1e-12,1e-3,100);
% loglog(xvec,xvec,'k--')
% loglog(xvec,xvec/10,'r--')
% loglog(xvec,xvec*10,'r--')
% 
% xlim([1e-12 1e-3])
% ylim([1e-12 1e-3])
% xlabel('log_{10}\chi [chameleon]','fontsize',16)
% ylabel('log_{10}\chi [CTD-\chi pod] ','fontsize',16)
% title('EQ14 CTD-\chi pod vs. chameleon mean')
% 
% if saveplot==1
%     SetPaperFigPath
%     figname=['EQ14_ctdChipod_vs_chamMean_chi_scatter']
%     print( fullfile( figdir , figname ), '-dpng' )
% end
%%