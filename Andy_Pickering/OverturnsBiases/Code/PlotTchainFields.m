%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotTchainFields.m
%
%
% 26 Nov. 2014 - AP - apickering@apl.washington.edu
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot temperature

clear ; close all

saveplot=1
plotiso=3

whmoor=4 % mooring #

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps'] ) )

clev=[1:2:8]
clev=[1:8]

Np=length(xx2.yday)

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,4);

cl=[1 8]

id=1:Np/4;

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

axes(ax(1))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
colorbar
caxis(cl)
title(['T-chain ' num2str(whmoor) ],'fontsize',16)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

id=round(2*Np/4) : 3*Np/4;

axes(ax(2))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
colorbar
caxis(cl)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

%
id=round(2*Np/4) : 3*Np/4;

axes(ax(3))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
colorbar
caxis(cl)
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

id=round(3*Np/4) : Np ;


axes(ax(4))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id))
cb=colorbar
ylabel(cb,'^oC')
caxis(cl)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)

linkaxes(ax,'y')
%

if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_TempOverview'])
    %    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    %    export_fig(fname,'-pdf','-r200')
    
    rmpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    save2pdfAP(fname)
    
end


%% Plot epsilon

clear ; close all
saveplot=1
plotiso=1

whmoor=1 % mooring #

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps']) )

Np=length(xx2.yday)

clev=[1:2:8]
clev=[1:8]


figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,4);

cl=[-9 -4]

id=1:Np/4;

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

axes(ax(1))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
colorbar
caxis(cl)
title(['T-chain ' num2str(whmoor) ],'fontsize',16)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

id=round(2*Np/4) : 3*Np/4;

axes(ax(2))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
colorbar
caxis(cl)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

%
id=round(2*Np/4) : 3*Np/4;

axes(ax(3))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
colorbar
caxis(cl)
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

id=round(3*Np/4) : Np ;


axes(ax(4))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar
ylabel(cb,'^log_{10}\epsilon')
caxis(cl)

if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),clev,'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)

linkaxes(ax,'y')
%

if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_EpsilonOverview'])
    %    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    %    export_fig(fname,'-pdf','-r200')
    
    rmpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    save2pdfAP(fname)
    
end
%

%%