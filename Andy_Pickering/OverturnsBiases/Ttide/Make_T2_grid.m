%~~~~~~~~~~~~~~~~~~~~~~~
%
% Make_T2_grid.m
%
% Interpolate T2 to evenly spaced depth grid and see how this affects
% overturns etc.
%
% See also ExamineT2.m
%
% 2 Mar 2015 - AP
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

dzgrid=10  % depth spacing in meters
dtgrid=10  % time spacing in minutes

load('/Users/Andy/Dropbox/ttide_jn/T2_short_timeseries.mat')

id=find(dat.T>3.5);
dat.T(id)=nan;

[Nz,Nt]=size(dat.T)

zgrid=dat.depth(1):-dzgrid:dat.depth(end);

Ti=nan*ones(length(zgrid),Nt);

for wht=1:Nt
    Ti(:,wht)=interp1(dat.depth,dat.T(:,wht),zgrid);
end

%%


T2=struct();
T2.yday=datenum2yday(dat.time);;
T2.z=1900-zgrid;
T2.t=Ti;
T2.Source='/Users/Andy/Dropbox/ttide_jn/T2_short_timeseries.mat'
T2.MakeInfo=['Made ' datestr(now) ' by AP w/ Make_T2_grid.m, in Matlab ' version]
T2.Info='Interpolated to evenly spaced depth grid, original time'
T2.dzgrid=dzgrid;
%T2.dtgrid=dtgrid;

Bdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide'
fname=fullfile(Bdir,'Data','T2',['T2_Gridded_dz_' num2str(dzgrid) 'm_traw'])
T2.fname=fname;
save(fname,'T2')

%%

figure(1);clf
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,2);

axes(ax(1))
ezpc(datenum2yday(dat.time),dat.depth,dat.T)
hold on
plot(datenum2yday(dat.time(1)),dat.depth,'ko')
colorbar
caxis([ 2.2128    3.4949])
axis xy

axes(ax(2))
ezpc(datenum2yday(dat.time),zgrid,Ti)
hold on
plot(datenum2yday(dat.time(1)),zgrid,'ko')
colorbar
caxis([ 2.2128    3.4949])
axis xy

linkaxes(ax)

%% interpolate to a coarser time grid also

yday0=datenum2yday(dat.time);
tgrid=nanmin(yday0) : dtgrid/60/24 : nanmax(yday0) ;
Ti2=nan*ones(length(zgrid),length(tgrid));
for whz=1:length(zgrid)
   Ti2(whz,:)=interp1(yday0,Ti(whz,:),tgrid);
end

%%

figure(1);clf
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,2);

axes(ax(1))
ezpc(datenum2yday(dat.time),zgrid,Ti)
hold on
%plot(datenum2yday(dat.time(1)),zgrid,'ko')
colorbar
caxis([ 2.2128    3.4949])
axis xy

axes(ax(2))
ezpc(tgrid,zgrid,Ti2)
hold on
%plot(datenum2yday(dat.time(1)),zgrid,'ko')
colorbar
caxis([ 2.2128    3.4949])
axis xy

linkaxes(ax)

%% save new gridded structure

T2=struct();
T2.yday=tgrid;
T2.z=1900-zgrid;
T2.t=Ti2;
T2.Source='/Users/Andy/Dropbox/ttide_jn/T2_short_timeseries.mat'
T2.MakeInfo=['Made ' datestr(now) ' by AP w/ Make_T2_grid.m, in Matlab ' version]
T2.Info='Interpolated to evenly spaced depth and time grid'
T2.dzgrid=dzgrid;
T2.dtgrid=dtgrid;

Bdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide'
fname=fullfile(Bdir,'Data','T2',['T2_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min'])
T2.fname=fname;
save(fname,'T2')

%%