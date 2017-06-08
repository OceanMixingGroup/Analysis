%
%
% Misc10Nov.m
%
% Similar to Misc8Nov.m but for up profiles only and sample every profile
%
% Compare with JN's results
%
%
% 10 Nov.
%
%%

clear ; close all
%whcase=1
tstart=tic

% now try resampling a densly sampled (in time) data set

disp('loading data to resample')
%~~ use JN's T-chains during IWISE
load('all_moorings.mat')
load('all_gridded.mat')
c=[4];
xx=grd{c}
clear grd


tshift=(2/60/24) % shift profiles in time by this much
Nprof=5040;%720

t_start=178%178%:tshift:177+Nshifts*tshift
close all
clear t_range z_range w_samp t_grd t_samp z_samp t0 z0 data_resamp

%    t_range=[t_start t_start+2] % ~time range to sample (days)
z_range=[580 2080] % depth-range to sample (m)
z_range=[xx.z(1) xx.z(end)]
%~~ vertical velocity of profiler
w_samp=0.25 % m/s

%~~ starting depth and time of profiler
%t0=t_range(1)
t0=t_start
z0=z_range(2);

dt_samp=40 % sampling dt (sec)
dz_samp=w_samp*dt_samp % sampling dz (m)

%~~~~~~ Everything below calculated from above chosen parameters

% time to do 1 profile
T_samp=range(z_range)/w_samp; % sec
T_samp_mins=T_samp/60; % min
T_samp_day=T_samp/86400; %days


%     % time of 'gridded MP' profiles==center of each profile
%     % each angled (t-z) profile will get gridded into one 'vertical' profile at
%     % the center time
%     t_grid0=t0+T_samp_day/2
%     t_grid=t_grid0  : T_samp_day : t_grid0+Nprof*T_samp_day    ;

%     % begin and end times of each profile
%     t_beg=nan*ones(Nprof,1);
%     t_end=nan*ones(Nprof,1);
%     t_beg=t0: T_samp_day : t0 + Nprof *T_samp_day;
%     t_end=t0+T_samp_day :T_samp_day : t0 + Nprof*T_samp_day +T_samp_day;

z_samp=[];

nt_per_prof=T_samp/dt_samp;
nt=nt_per_prof*Nprof;

% new calc of depth-time path of profiler
% depth vector is same for each, so only make once
% also only make tvec once, then shift for each profile
%~~~

z_samp = z_range(1) : dz_samp : z_range(2);
z_samp=flipud(z_samp(:));

dt_samp_day=dt_samp/86400;

t_samp1=t0 : dt_samp_day :  t0 + dt_samp_day*length(z_samp);
%

figure(1);clf
plot(t_samp1(1:length(z_samp)),z_samp,'o-')
hold on
axis ij
ylabel('Depth (m) ')
ylabel('Time (yearday)')
%

%
t_real=datenum2yday(xx.time);
z_real=xx.z;
data_real=xx.T;

data_resamp=nan*ones(length(z_samp),Nprof);
t_grid=nan*ones(1,Nprof);
hb=waitbar(0,'resampling data')
for whp=1:Nprof
    waitbar(whp/Nprof,hb)
    clear t_samp
    t_samp=t_samp1+ whp*tshift;
    
    t_grid(whp)=nanmean(t_samp);
    % crudely  'resample' by finding real data points closest to each sampling point in
    % depth/time.  In future, maybe use interp2?
    %    disp('resampling data')
    
    for whz=1:length(z_samp)
        
        clear Iz It valz
        [valz,Iz]=nanmin(abs(z_samp(whz)-z_real));
        [valt,It]=nanmin(abs(t_samp(whz)-t_real));
        
        if valz<10 && valt<2/60/24
            data_resamp(whz,whp)=data_real(Iz,It);
            %        figure(1)
            %      plot(t_real(It),z_real(Iz),'x')
        else
            disp('adfl')
            
        end
    end
    
    
end
delete(hb)

%%
%
idtr=isin(t_real,[t_grid(1) t_grid(end)]);

xl=[t_grid(1) t_grid(end)]

figure(1);clf

cl=[1 8]
yl=z_range
figure(3);clf
agutwocolumn(0.8)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,2);

axes(ax(1));
ezpc(t_real(idtr),z_real,data_real(:,idtr))
hold on
contour(t_real(idtr),z_real,data_real(:,idtr),[0:0.5:6],'k')
%    plot(t_samp,z_samp,'w-','linewidth',1)
ylim(yl)
xlim(xl)
caxis(cl)
colorbar
SubplotLetterMW('Actual')

axes(ax(2))
ezpc(t_grid,z_samp,data_resamp)
hold on
contour(t_grid,z_samp,data_resamp,[0:0.5:6],'k')
ylim(yl)
xlim(xl)
caxis(cl)
colorbar
SubplotLetterMW('Resampled')

linkaxes(ax)

%% now compute overturns for the resampled profiles


%~~ T-chain dz is 10 m
minOT=1; % min size of overturn

%~~ First do resampled data

x_resamp.time=t_grid;
x_resamp.z=flipud(z_samp);
x_resamp.T=flipud(data_resamp);
x_resamp.S=34.604-.045*(x_resamp.T-2.5);
x_resamp.eps=NaN*x_resamp.S;
x_resamp.Lot=NaN*x_resamp.S;
x_resamp.Lttot=NaN*x_resamp.S;
%
addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/

h = waitbar(0,'Please wait...');
%
for a=1:length(x_resamp.time)
    ind=a;
    if mean(x_resamp.T(:,ind)<6)
        %        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.P(:,ind),xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(x_resamp.z(:),x_resamp.T(:,ind),x_resamp.S(:,ind),35.8,0,minOT,1e-5,0);
        x_resamp.eps(:,ind)=Epsout;
        x_resamp.Lot(:,ind)=Lot;
        x_resamp.Lttot(:,ind)=Lttot;
    end
    %    if mod(a,1000)==0
    waitbar(a/length(x_resamp.time),h)
    %  end
end

delete(h)

%% also recompute overturns for real data (there is already an eps field, but get different answer when recomutpe...)

xx.yday=datenum2yday(xx.time);
xx2=xx;
xx2.time=xx.time(idtr);
xx2.yday=xx.yday(idtr);
xx2.T=xx.T(:,idtr);
xx2.S=34.604-.045*(xx2.T-2.5);
xx2.eps=NaN*xx2.S;
xx2.Lot=NaN*xx2.S;
xx2.Lttot=NaN*xx2.S;

hb=waitbar(0)
for ind=1:length(xx2.time)
    waitbar(ind/length(xx2.time),hb)
        if mean(xx2.T(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,1,1e-5,0);
        xx2.eps(:,ind)=Epsout;
        xx2.Lot(:,ind)=Lot;
                xx2.Lttot(:,ind)=Lttot;
    end
  
end

delete(hb)
%%
figure(2);clf
agutwocolumn(0.8)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,2);

axes(ax(1))
%ezpc(t_real(idtr),z_real,log10(xx.eps(:,idtr)))
ezpc(xx2.yday,xx2.z,log10(xx2.eps))
hold on
contour(t_real(idtr),z_real,data_real(:,idtr),[0:0.5:6],'k')
caxis([-8 -3])
colorbar

axes(ax(2));
ezpc(x_resamp.time,x_resamp.z,log10(x_resamp.eps))
hold on
contour(t_grid,z_samp,data_resamp,[0:0.5:6],'k')
caxis([-8 -3])
colorbar
cmap=jet;
cmap=[0.7*[1 1 1];cmap;]
colormap(cmap)
SubplotLetterMW('resamp')
linkaxes(ax)

%%

figure(3);clf

semilogx(nanmean(x_resamp.eps,2),x_resamp.z,'k','linewidth',2)
hold on
%semilogx(nanmean(xx.eps(:,idtr),2),z_real,'m','linewidth',2)
semilogx(nanmean(xx2.eps,2),xx2.z,'r','linewidth',2)
axis ij


%%

figure(22);clf
subplot(121)
hist(xx2.Lot(:),30)

subplot(122)
hist(x_resamp.Lot(:),30)


%%



figure(22);clf
subplot(121)
hist(xx2.Lttot(:),30)
xlim([0 300])
subplot(122)
hist(x_resamp.Lttot(:),30)
xlim([0 300])
%%

figure(22);clf
subplot(121)
hist(log10(xx2.eps(:)),30)
xlim([-9 -4])

subplot(122)
hist(log10(x_resamp.eps(:)),30)
xlim([-9 -4])

%%