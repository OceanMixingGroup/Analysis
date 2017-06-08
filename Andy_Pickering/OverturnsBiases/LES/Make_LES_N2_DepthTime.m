%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_LES_N2_DepthTime.m
%
% Read data from LES model at N2 location, run by Sutanu and
% student Masoud, ad make into a gridded depth-time structure.
%
% Formerly part of PlotLES_N2.m
%
% 20 Feb 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_all.mat')

% Make a depth-time matrix to pcolor

indy=30;

allu=nan*ones(193,1759);
allth=nan*ones(193,1759);
allw=nan*ones(193,1759);
yday=nan(ones(1,1759));
for ind=1:1759
    allu(:,ind)=interlav(ind).u(indy,:);
    allv(:,ind)=interlav(ind).v(indy,:);
    allw(:,ind)=interlav(ind).w(indy,:);
    allth(:,ind)=interlav(ind).th(indy,:);
    yday(ind)=interlav(ind).year_time;
end

N2=struct();
N2.Yvalue=indy;
N2.u=flipud(allu);
N2.v=flipud(allv);
N2.w=flipud(allw);
N2.sgth=flipud(allth*100 +1026);
N2.yday=yday;
z=abs(interlav(1).z(1,:))*100;z=z(:);
N2.z=flipud(z);
N2.MakeInfo=['Made ' datestr(now) ' w/ Make_LES_N2_DepthTime.m ']
fname='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_depth_time';
N2.fname=fname
save(fname,'N2')

%% Plot summary of data

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,4);

thm=nanmean(N2.sgth,2);
zc=0:10:nanmax(N2.z);
ti=interp1(N2.z,thm,zc);

dc=10

axes(ax(1))
ezpc(N2.yday,N2.z,N2.u)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
%axis xy
colorbar
caxis(0.5*[-1 1])
SubplotLetterMW('u')
ylabel('Depth','fontsize',16)
%xtloff

axes(ax(2))
ezpc(N2.yday,N2.z,N2.v)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
colorbar
caxis(0.1*[-1 1])
SubplotLetterMW('v')
ylabel('Depth','fontsize',16)

axes(ax(3))
ezpc(N2.yday,N2.z,N2.w)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
colorbar
SubplotLetterMW('w')
ylabel('Depth','fontsize',16)
caxis(0.1*[-1 1])
colormap(bluered)
%xtloff


axes(ax(4))
ezpc(N2.yday,N2.z,N2.sgth)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
colorbar
%caxis(0.5*[-1 1])
SubplotLetterMW('\rho')
%xtloff
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
colormap(ax(4),jet)

linkaxes(ax)

%
%% Also make a structure with data interpolated to evenly spaced depth and time grid

zc=0:5:nanmax(N2.z);

[M,N]=size(N2.u);

ui=nan*ones(length(zc),length(N2.yday));
vi=nan*ones(length(zc),length(N2.yday));
wi=nan*ones(length(zc),length(N2.yday));
sgthi=nan*ones(length(zc),length(N2.yday));

for wht=1:N
    ui(:,wht)=interp1(N2.z,N2.u(:,wht),zc);
    vi(:,wht)=interp1(N2.z,N2.v(:,wht),zc);
    wi(:,wht)=interp1(N2.z,N2.w(:,wht),zc);
    sgthi(:,wht)=interp1(N2.z,N2.sgth(:,wht),zc);
end

N2even=struct();
N2even.u=ui;
N2even.v=vi;
N2even.w=wi;
N2even.sgth=sgthi;
N2even.z=zc;
N2even.yday=N2.yday;

N2even.MakeInfo=['Made ' datestr(now) ' w/ Make_LES_N2_DepthTime.m ']
fname='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_depth_time_Even';
N2even.fname=fname
save(fname,'N2even')

%% Make a summary figure

clear ; close all

saveplot=1

load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_depth_time_Even')

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,4);

thm=nanmean(N2even.sgth,2);

dc=20;

axes(ax(1))
ezpc(N2even.yday,N2even.z,N2even.u)
hold on
contour(N2even.yday,N2even.z,N2even.sgth,thm(1:dc:end),'k')
colorbar
caxis(0.5*[-1 1])
SubplotLetterMW('u')
ylabel('Depth','fontsize',16)
title('LES Model at N2')

axes(ax(2))
ezpc(N2even.yday,N2even.z,N2even.v)
hold on
contour(N2even.yday,N2even.z,N2even.sgth,thm(1:dc:end),'k')
colorbar
caxis(0.1*[-1 1])
SubplotLetterMW('v')
ylabel('Depth','fontsize',16)

axes(ax(3))
ezpc(N2even.yday,N2even.z,N2even.w)
hold on
contour(N2even.yday,N2even.z,N2even.sgth,thm(1:dc:end),'k')
colorbar
SubplotLetterMW('w')
ylabel('Depth','fontsize',16)
caxis(0.1*[-1 1])
colormap(bluered)

axes(ax(4))
ezpc(N2even.yday,N2even.z,N2even.sgth)
hold on
contour(N2even.yday,N2even.z,N2even.sgth,thm(1:dc:end),'k')
colorbar
SubplotLetterMW('\rho')
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
colormap(ax(4),tempmap);

linkaxes(ax)

if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['LES_N2_summary'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-png')
end

%%