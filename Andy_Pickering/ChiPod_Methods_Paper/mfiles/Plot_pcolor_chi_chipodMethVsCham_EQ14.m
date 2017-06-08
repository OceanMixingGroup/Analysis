%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_pcolor_chi_chipodMethVsCham_EQ14.m
%
% * Makes plot for chipod methods paper *
%
% Was part of Plot_chiProc_vs_actual_EQ14.m before
%
%---------------
% 04/06/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

saveplot=1

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

% Set chipod params and load chameleon data
SetParamsPaperPlots
Params

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/',...
    ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]));

% average chameleon data in same depth bins as chipod method
% so we can make scatter plots etc further on
cham_bin=MakeBinnedChamEq14(C.BinParams,cham)

% Pcolor chi comparison

xl=[370 3100];
yl=[0 200];
cl=[-10 -5];

id1=find(C.castnumber==374)
id2=find(cham.castnumber==374)

figure(2);clf
agutwocolumn(0.8)
wysiwyg

ax1=subplot(211);
ezpc(C.castnumber(id1:end),sw_dpth(C.p,0),log10(C.chi(:,id1:end)))
cb=colorbar;
cb.Label.String='log_{10}\chi [K^2s^{-1}]';cb.Label.FontSize=14;
xlim(xl)
ylim(yl)
caxis(cl)
title('\chi pod method \chi_{\chi}^{cham}','fontsize',16);
ylabel('Depth [m]','fontsize',16)
xlabel('Cast #','fontsize',16)

chamnames = ['0004';'0400';'0453';'0507';'0588';'0642';'0705';'0753';...
    '0754';'0794';'0873';'0957';'1040';'1135';'1301';'1405';...
    '1533';'1590';'1591';'1652';'1795';'1903';'2003';'2099';...
    '2196';'2301';'2390';'2479';'2581';'2670';'2761';'2853';...
    '2952';'3029';'3089'];

hold on
for i=1:length(chamnames)
    chamstr=chamnames(i,:);
    cnum1=str2num(chamstr);
    plot(cnum1,0,'kd')
end

ax2=subplot(212);
ezpc(cham.castnumber(id2:end),sw_dpth(cham.P(:,id2:end),0),log10(cham.CHI(:,id2:end)))
cb=colorbar;
cb.Label.String='log_{10}\chi [K^2s^{-1}]';cb.Label.FontSize=14;
xlim(xl)
ylim(yl)
caxis(cl)
title('Chameleon \chi_{\epsilon}^{cham}','fontsize',16);
xlabel('Cast #','fontsize',16)
ylabel('Depth [m]','fontsize',16)

linkaxes([ax1 ax2])

if saveplot==1
    figname=['EQ14_chiVsTrue_chi_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
    SetPaperFigPath
    print( fullfile( figdir,figname) , '-dpng' )
end

%%