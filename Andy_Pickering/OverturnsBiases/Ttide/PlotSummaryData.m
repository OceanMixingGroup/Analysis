
%%

clear ; close all

minOT=10

DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT)]

load(fullfile(DatDir,fname))

%%

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,2);

tm=nanmean(xx2.T,2);
zc=xx2.z(1):10:xx2.z(end);
tmi=interp1(xx2.z,tm,zc)

axes(ax(1))
ezpc(xx2.yday,xx2.z,xx2.T)
hold on
contour(xx2.yday,xx2.z,xx2.T,tmi(1:10:end),'k')
cb=colorbar
cb.Label.String='^oC';
colormap(ax(1),'parula')
title('T-tide Mooring T2 ')

axes(ax(2))
ezpc(xx2.yday,xx2.z,log10(xx2.eps))
hold on
contour(xx2.yday,xx2.z,xx2.T,tmi(1:10:end),'k')
cb=colorbar
cb.Label.String='log_{10}\epsilon';
caxis([-9 -5])
cmap=jet;
colormap(ax(2),[ 0.8*[1 1 1] ; cmap])
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
linkaxes(ax)


fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['T2_Temp_and_Eps'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-png')

%%

h=OT_Hist_Summary(xx2)

fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['T2_OThist_raw'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-png')
