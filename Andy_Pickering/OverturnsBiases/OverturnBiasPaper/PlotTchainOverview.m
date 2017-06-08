%
%
%
% June 17, 2015 - A. Pickering
%
%% Plot epsilon

clear ; close all
saveplot=1
plotiso=1
minOT=50
whmoor=3 % mooring #

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data',['Tchain' num2str(whmoor)] , ['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]) )
%
%clev=[1:2:8]
clev=[0.5:8.5]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));


figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot2(0.1, 0.03, 0.02, 0.06, 0.1, 0.005, 1,3);

cl=[-9 -4]

xl=[175 190]
id=isin(xx2.yday,xl);

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

load('/Users/Andy/Cruises_Research/IWISE/Data/TPXO/TPXO_Iwise2011_IWISE11_N2.mat')

% plot BT velocity
axes(ax(1))
plot(TS.yday,TS.ALL.u/100)
xlim(xl)
gridxy
cb=colorbar;
killcolorbar(cb)
title(['T-chain ' num2str(whmoor) ],'fontsize',16)
ylabel('m/s','fontsize',16)
SubplotLetterMW('U_{BT}')
xtloff

axes(ax(2))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
colorbar
%caxis(cl)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
%       contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:50:end),'k')
end
SubplotLetterMW('T')
ylabel('Depth','fontsize',16)
xtloff


axes(ax(3))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
colorbar
caxis(cl)

if plotiso==1
    hold on
%    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
%       contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:50:end),'k')
end


cmap=flipud(hot);
colormap(gca,[0.75*[1 1 1] ; cmap])

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
SubplotLetterMW('\epsilon')

linkaxes(ax,'x')
%
if saveplot==1
%    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',)
    
        figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_Overview']
    print('-dpng',fullfile(figdir,figname))

end

%%