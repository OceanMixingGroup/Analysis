%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotHistChiRatio_CTDchamPairs.m
%
% * Makes plot for chipod methods paper *
%
% Plot histograms of the ratio of chipod to chameleon chi, as well as the
% ratio of before/after chameleon pairs, to show that difference between
% methods is similar to natural variability.
%
% See Make_Combined_Cham_for_CTD_pairs.m
%
%-------------------
% 05/12/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=1

load('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Ctd_Cham_Pairs.mat')


%% Plot hist of ratio of chis for cham before/after, and chipod-cham
% want to show that spread in chameleon profiles is similar

chamrat=Pairs.CHAM.after.chi(:) ./ Pairs.CHAM.before.chi(:);
chirat=Pairs.CTD.chi1(:) ./ Pairs.CHAM.before.chi(:) ;
chirat2=Pairs.CTD.chi1(:) ./ Pairs.CHAM.after.chi(:) ;


Ds='stair'
Ds='bar'
Nm='pdf'

figure(1);clf
agutwocolumn(0.6)
wysiwyg
set(gcf,'defaultaxesfontsize',14)
h1=histogram(log10(chamrat),'DisplayStyle',Ds,'FaceColor','r','Normalization',Nm)
hold on
h2=histogram(log10(chirat),h1.BinEdges,'DisplayStyle',Ds,'FaceColor','b','Normalization',Nm)
h3=histogram(log10(chirat2),h1.BinEdges,'DisplayStyle',Ds,'FaceColor','g','Normalization',Nm)
grid on
xlabel('log_{10}[\chi_1 / \chi_2]','fontsize',16)
hleg=legend([h1 h2 h3],'\chi_{\epsilon 1}^{cham} / \chi_{\epsilon 2}^{cham}','\chi_{\chi}^{ctd}/\chi_{\epsilon 1}^{cham}','\chi_{\chi}^{ctd}/\chi_{\epsilon 2}^{cham}')
hleg.FontSize=16;
ylabel('Pdf','fontsize',16)

freqline(nanmean(log10(chamrat)),'r--')
freqline(nanmean(log10(chirat)),'b--')
freqline(nanmean(log10(chirat2)),'g--')

xlim(4.2*[-1 1])

% Add some stats to figure

text(1.4,.3,['\mu = ' num2str(roundx(nanmean(log10(chamrat)),2)) ', \sigma= ' num2str(roundx(nanstd(log10(chamrat)),2))],'fontsize',17,'color','red')
text(1.4,.25,['\mu = ' num2str(roundx(nanmean(log10(chirat)),2)) ', \sigma= ' num2str(roundx(nanstd(log10(chirat)),2))],'fontsize',17,'color','b')
text(1.4,.2,['\mu = ' num2str(roundx(nanmean(log10(chirat2)),2)) ', \sigma= ' num2str(roundx(nanstd(log10(chirat2)),2))],'fontsize',17,'color',[0 0.5 0])%,'fontweight','bold')


%%

if saveplot==1
    
    SetPaperFigPath
    figname='EQ14_CtdChipod_hist_chi_ratios'
    print( fullfile( figdir , figname ) ,'-dpng')
    
end

%%