%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareLES_N2_OTvsDiss.m
%
%
% 24 Feb. 2015
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('N2diss.mat')
load('N2_OT.mat')

%%

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.075, 0.06, 0.1, 0.05, 1,2);

axes(ax(1))
ezpc(OT.yday,OT.z,log10(OT.eps))
%axis xy
caxis([-10 -4])
colorbar
cmap=jet;
colormap([ 0.9*[1 1 1] ; cmap])
ylabel('Depth')
SubplotLetterMW('OT')

axes(ax(2))
ezpc(N2diss.yday,N2diss.z,log10(N2diss.eps))
%ezpc(N2diss.yday,N2diss.z,(N2diss.eps))
%axis xy
caxis([-10 -4])
colorbar
xlabel('Yearday')
ylabel('Depth')
SubplotLetterMW('diss')

linkaxes(ax)

%%

saveplot=1

figure(2);clf
agutwocolumn(0.5)
wysiwyg
semilogx(nanmean(OT.eps,2),OT.z)
hold on
semilogx(nanmean(N2diss.eps,2),N2diss.z)
axis ij
legend('OT','Diss')
shg
grid on
xlim([1e-9 1e-4])
xlabel('<\epsilon> [Wkg^{-1}]')
ylabel('Depth [m]')
title('LES at N2')


if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['LES_N2_OTvsDirectProfiles'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-png')
end
%%
