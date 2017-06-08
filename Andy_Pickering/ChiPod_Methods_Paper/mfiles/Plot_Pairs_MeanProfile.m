%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Pairs_MeanProfile.m
%
% * Makes plot for chipod methods paper *
%
%
%-------------
% 05/12/15 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot time-mean profiles w/ bootstrap confidence limits

clear ; close all

saveplot=1

load('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Ctd_Cham_Pairs.mat')

% depth interval to average over
interval=50;

figure(1);clf
agutwocolumn(0.75)
wysiwyg

%[~, ~ , hcham]=BinAndBootProfiles([Pairs.CHAM.before.chi Pairs.CHAM.after.chi],Pairs.CHAM.p,interval,100,1,0);
[~, ~ , hcham]=BinAndBootProfiles([Pairs.CHAM.before.chi Pairs.CHAM.after.chi],sw_dpth(Pairs.CHAM.p,0),interval,100,1,0);
hold on
%semilogx(nanmean(Pairs.CTD.chi1,2),Pairs.CTD.p)
%[~ , ~ , hchi]=BinAndBootProfiles(Pairs.CTD.chi1,Pairs.CTD.p,interval,100,1,0);
[~ , ~ , hchi]=BinAndBootProfiles(Pairs.CTD.chi1,sw_dpth(Pairs.CTD.p,0),interval,100,1,0);
%semilogx(nanmean(cham_all.chi2,2),cham_all.p)
axis ij
grid on

ylabel('Depth [m]','fontsize',16)
xlabel('log_{10}\chi [K^2s^{-1}]','fontsize',16)
title(['EQ14 CTD-\chi pod vs Cham - All Pairs'])
legend([hcham hchi],'\chi_{\epsilon}^{cham}','\chi_{\chi}^{ctd}','location','best')

ylim([0 198])
xlim([1e-8 1e-5])
set(gca,'fontsize',14)
%ylabel('pressure [db]','fontsize',16)
pos=get(gca,'position');
set(gca,'position',pos.*[1 1.3 0.9 0.9])


if saveplot==1
    figname=['EQ14_chi_cham_meanProf_all']
    SetPaperFigPath
    print( fullfile( figdir , figname ) , '-dpng' )
end

%%