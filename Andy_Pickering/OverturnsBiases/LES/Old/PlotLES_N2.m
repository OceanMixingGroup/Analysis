%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotLES_N2.m
%
% Examine and plot data from LES model at N2 location, run by Sutanu and
% student Masoud.
%
%
% 15 Feb 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_all.mat')

%% Plot some profiles in time

for ind=1:5:500
    figure(1);clf
    plot(interlav(ind).th,interlav(ind).z(1,:)*100)
    %xlim(0.2*[-1 1])
    pause(1)
end

%% Make a depth-time matrix to pcolor

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
N2.u=allu;
N2.v=allv;
N2.w=allw;
N2.sgth=allth*100 +1026;
N2.yday=yday;
N2.z=interlav(1).z(1,:)*100;

fname='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_depth_time';
save(fname,'N2')

%%

close all

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,4);

thm=nanmean(N2.sgth,2);
zc=nanmin(N2.z):10:0;
ti=interp1(N2.z,thm,zc);

dc=10

axes(ax(1))
ezpc(N2.yday,N2.z,N2.u)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
axis xy
colorbar
caxis(0.5*[-1 1])
SubplotLetterMW('u')
ylabel('Depth','fontsize',16)
%xtloff

axes(ax(2))
ezpc(N2.yday,N2.z,N2.v)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
axis xy
colorbar
caxis(0.1*[-1 1])
SubplotLetterMW('v')
%xtloff
ylabel('Depth','fontsize',16)

axes(ax(3))
ezpc(N2.yday,N2.z,N2.w)
hold on
contour(N2.yday,N2.z,N2.sgth,ti(1:dc:end),'k')
axis xy
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
axis xy
colorbar
%caxis(0.5*[-1 1])
SubplotLetterMW('\rho')
%xtloff
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
colormap(ax(4),jet)

linkaxes(ax)

%%