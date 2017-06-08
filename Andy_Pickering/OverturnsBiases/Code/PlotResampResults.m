%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotResampResults.m
%
% Plot results of resampling T-chains (done in ResampleProfiler_V2.m)
%
%
%
% Originally modified from PlotResampResults12Nov.m
%
% 24 Feb. 2015 - A. Pickering - andypicke@gmail.com
% 28 May 2015  - Update for new format
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot mean depth-profiles of epsilon from 1 resampling set (1 speed)

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3;xl=[1e-9 10^(-5.5)]
%whmoor=4;xl=[1e-9 10^(-6.5)]
%whmoor=1;xl=[1e-9 10^(-7)]
%whmoor=2;xl=[1e-9 10^(-6)]

minOT=50
testnum=1 % 0.25m/s
%testnum=2 % 0.5m/s
% testnum=3 % 0.15m/s
%testnum=4 % 0.75m/s
testnum=6


figure(1);clf
ezall=[];

% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%~~  Load 'true' T-chain data 
load(fullfile('Data',['Tchain' num2str(whmoor)],['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

% find data within time of resampling period
%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);

Nshift=size(REsamp.eps,3)
ezall=nan*ones(length(REsamp.z),Nshift);
for whcase=1:Nshift
    ezall(:,whcase)=nanmean(REsamp.eps(:,:,whcase),2);
   semilogx(nanmean(REsamp.eps(:,:,whcase),2),REsamp.z,'.','color',0.7*[1 1 1])
   hold on
end
axis ij

% plot actual time-mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)

% plot resampled time-mean profile
hs=semilogx(nanmean(ezall,2),REsamp.z,'linewidth',3)

axis ij
xlim(xl)
ylim([xx2.z(1) xx2.z(end)])
xlabel('<\epsilon> [Wkg^{-1}]','fontsize',16)
ylabel('Depth [m]','fontsize',16)
title(['Tchain ' num2str(whmoor) ' w=' num2str(REsamp.w_samp) ' ms^{-1}'])
legend([ht hs],'True','Sampled')
grid on
shg
%
%
% figdir='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/'
% fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth' ])
% save2pdfAP(fname)

%% same as above, but bin average and bootstrap

figure(2);clf
interval=100
nboot=100
addpath /Users/Andy/Cruises_Research/ChiPod/IWISE10_And_vmp/
[xout bb hs]=BinAndBootProfiles(ezall,REsamp.z,interval,nboot,1);
hold on
[xout bb ht]=BinAndBootProfiles(xx2.eps(:,idtr),xx2.z,interval,nboot,1);
%ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
legend([hs ht])


%% Plot mean depth-profiles of epsilon for resampling at different sampling speeds 

clear ; close all

saveplot=0

whmoor=3
minOT=50

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases


%~~  Load T-chain data to resample
% load ' true ' datea
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
%load(fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);


ezall1=[];
testnum=1
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall1=squeeze(nanmean(REsamp.eps,2));

zvec=REsamp.z;
%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);

clear REsamp

ezall2=[];
testnum=2
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall2=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall3=[];
testnum=3
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall3=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall4=[];
testnum=4
% load resampled data
clear REsamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall4=squeeze(nanmean(REsamp.eps,2));
clear REsamp

figure(1);clf
agutwocolumn(0.6)
wysiwyg

% plot resampled mean profiles
h1=semilogx(nanmean(ezall1,2),zvec,'linewidth',2)
hold on
h2=semilogx(nanmean(ezall2,2),zvec,'linewidth',2)
h3=semilogx(nanmean(ezall3,2),zvec,'linewidth',2)
h4=semilogx(nanmean(ezall4,2),zvec,'linewidth',2)

% plot actual mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
axis ij
xlim([1e-9 10^(-5.5)])
xlabel('<\epsilon> [Wkg^{-1}]','fontsize',16)
ylabel('Depth [m]','fontsize',16)
title(['Tchain ' num2str(whmoor)])
grid on
shg

legend([ht h3 h1 h2 h4 ],'Actual','15cm/s','25cm/s','50cm/s','75cm/s','location','best')

%
if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth_DiffSpeeds' ])
    export_fig(fname,'-pdf','-r200')
end

%
%
%% Plot mean depth-profiles of epsilon for resampling at different sampling speeds 
% ** 2 X 2 figure, each panel is different w_samp
% ** plot all ensembles in addition to mean
%  25 Feb - update to use REsamp etc.

clear ; close all

saveplot=1
minOT=50
whmoor=3 ; xl=[1e-9 10^(-5.5)]
%whmoor=4 ; xl=[1e-9 10^(-5.9)]
%whmoor=1 ; xl=[1e-9 10^(-5.9)]
%whmoor=2 ; xl=[1e-9 10^(-5.9)]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
%

ezall1=[];
testnum=1
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall1=squeeze(nanmean(REsamp.eps,2));

%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);
zvec=REsamp.z;
clear REsamp

ezall2=[];
testnum=2
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall2=squeeze(nanmean(REsamp.eps,2));
clear REsamp

ezall3=[];
testnum=3
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall3=squeeze(nanmean(REsamp.eps,2));
clear REsamp

%
ezall4=[];
testnum=4
% load resampled data
clear REsamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall4=squeeze(nanmean(REsamp.eps,2));
clear REsamp

%
myc=[0.8500    0.3250    0.0980];
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 2,2);

axes(ax(1))
% plot resampled mean profiles
semilogx(ezall3,zvec,'.','color',0.7*[1 1 1])
hold on
h3=semilogx(nanmean(ezall3,2),zvec,'linewidth',2,'color',myc)
hold on
% plot actual mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% plot +/- standard devation
semilogx(nanmean(ezall3,2)-std(ezall3,0,2),zvec,'--','color',0.3*[1 1 1])
semilogx(nanmean(ezall3,2)+std(ezall3,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
axis ij
ylabel('Depth ','fontsize',16)
ht=SubplotLetter(['w=0.15 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on

axes(ax(2))
% plot resampled mean profiles
semilogx(ezall1,zvec,'.','color',0.7*[1 1 1])
hold on
h1=semilogx(nanmean(ezall1,2),zvec,'linewidth',2,'color',myc)
hold on
% plot actual mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% plot +/- standard devation
semilogx(nanmean(ezall1,2)-std(ezall1,0,2),zvec,'--','color',0.3*[1 1 1])
semilogx(nanmean(ezall1,2)+std(ezall1,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
ytloff
axis ij
ht=SubplotLetter(['w=0.25 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on

axes(ax(3))
semilogx(ezall2,zvec,'.','color',0.7*[1 1 1])
hold on
% plot resampled mean profiles
h2=semilogx(nanmean(ezall2,2),zvec,'linewidth',2,'color',myc)
hold on
% plot actual mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% plot +/- standard devation
semilogx(nanmean(ezall2,2)-std(ezall2,0,2),zvec,'--','color',0.3*[1 1 1])
semilogx(nanmean(ezall2,2)+std(ezall2,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
axis ij
xlabel('<\epsilon> (Wkg^{-1})','fontsize',16)
ylabel('Depth ','fontsize',16)
ht=SubplotLetter(['w=0.5 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on

axes(ax(4))
semilogx(ezall4,zvec,'.','color',0.7*[1 1 1])
hold on
% plot resampled mean profiles
h4=semilogx(nanmean(ezall4,2),zvec,'linewidth',2,'color',myc)
hold on
% plot actual mean profile
ht=semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',3)
% plot +/- standard devation
semilogx(nanmean(ezall4,2)-std(ezall4,0,2),zvec,'--','color',0.3*[1 1 1])
semilogx(nanmean(ezall4,2)+std(ezall4,0,2),zvec,'--','color',0.3*[1 1 1])
xlim(xl)
ylim([zvec(1) zvec(end)])
ytloff
axis ij
xlabel('<\epsilon> (Wkg^{-1})','fontsize',16)
ht=SubplotLetter(['w=0.75 m/s'],0.6,0.1)
set(ht,'fontsize',14)
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
grid on
shg

linkaxes(ax)

if saveplot==1    
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_Resamp_EpsvsDepth_DiffSpeeds_2X2' ]
    print('-dpng',fullfile(figdir,figname))    
end

%
%
%% Plot depth-integrated  epsilon for 1 resample set (1 speed)

clear ; close all

saveplot=0

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/

minOT=50
t_avg=24/24
whmoor=3
testnum=1

% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])

clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
ezall1=squeeze(nanmean(REsamp.eps,2));

%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);
zvec=REsamp.z;

eall_re=[];
tall_re=[];

figure(1) ; clf

Nshift=size(REsamp.eps,3)
for whcase=1:Nshift
    
    clear epsint_re x_resamp eps_re t_re
    
    dz_samp=nanmean(diff(REsamp.z));
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
%    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.timeall(whcase,:) ) ;
    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;
    semilogy(t_re,eps_re,'.','color',0.75*[1 1 1],'linewidth',1)
    hold on
    
    eall_re=[eall_re; eps_re];
    tall_re=[tall_re ; t_re];
    
end

%
epsint=nansum(xx2.eps)*10*1026;
[eps,t]=SimpleBoxCar(epsint,t_avg,xx2.yday);

semilogy(t,eps,'k','linewidth',3)

semilogy(t_re,nanmean(eall_re),'linewidth',2)%,'color',0.5*[1 1 1])

grid on
title(['Tchain ' num2str(whmoor) ' \int \epsilon dz, wsamp=' num2str(REsamp.w_samp) 'm/s , ' num2str(t_avg) ' day avg. '])
xlabel('Yearday','fontsize',16)
ylabel('\int \rho \epsilon dz / Wm^{-2} ','fontsize',16)

%xlim([165 185])
ylim([1e-4 1e1])

shg
%
if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    
    fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Resamp_IntEpsvsTime' ])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
    
end

%% make a scatter plot

% 22 Dec. 2014

saveplot=0

figure(3);clf

% t,eps  ,  t_re, ,nanmean(eall_re)
% interp to same times

ei2=interp1(t_re,nanmean(eall_re),t);
bias=ei2-eps;
%bias=(ei2-eps)./eps*100;
%semilogx(eps,100*bias./eps,'o','linewidth',2)
loglog(eps,abs(bias),'p','linewidth',2)
grid on
xlabel('\epsilon_{true}')
ylabel('\epsilon_{resamp}-\epsilon_{true}')
%xvec=linspace(nanmin(eps),nanmax(eps),50);
xvec=linspace(1e-6,1e2,50);
hold on
loglog(xvec,xvec,'k--')
title(['Tchain ' num2str(whmoor) ' - \int \epsilon dz - ' num2str(t_avg*24) ' hr avg' ])

SubplotLetter(['w_{samp}=' num2str(REsamp.w_samp) 'm/s'])

if saveplot==1
figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
fname=fullfile(figdir,['Tchain' num2str(whmoor) '_IntEps_ErrorvsTrue' ])
rmpath '/Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/'
save2pdfAP(fname)
end
%

figure(33);clf
subplot(211)
ezpc(xx2.yday,xx2.z,xx2.T)
colorbar
subplot(212)
ezpc(xx2.yday,xx2.z,log10(xx2.eps))
colorbar
caxis([-9 -4])


%
%% Plot depth-integrated epsilon timeseries for 4 different sampling speeds
% ** make a 4 panel plot, each panel is for a different profilig speed

% **** 25 Feb - need to interp to common time vector to average ensembles
%

clear ; close all

saveplot=1

whmoor=3 ; yl=[10^(-3.5) 1e1] ;
%whmoor=4 ; yl=[ 10^(-3.5) 10^(0.5)] ;
%whmoor=1 ; yl=[ 10^(-4.5) 10^(0.5)] ;
t_avg=1
minOT=50

%~~
testnum=1
clear eall_re tall_re
eall_re=[];
tall_re=[];
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Nshift=size(REsamp.eps,3)
for whcase=1:Nshift    
    clear epsint_re x_resamp eps_re t_re    
    dz_samp=nanmean(diff(REsamp.z));
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
%    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.timeall(whcase,:) ) ;
     [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;
    eall_re=[eall_re; eps_re];
    tall_re=[tall_re ; t_re];
end
eall_re1=eall_re;
tall_re1=tall_re;

%~~
testnum=2
clear eall_re tall_re
eall_re=[];
tall_re=[];
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Nshift=size(REsamp.eps,3)
for whcase=1:Nshift    
    clear epsint_re x_resamp eps_re t_re    
    dz_samp=nanmean(diff(REsamp.z));
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
%    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.timeall(whcase,:) ) ;
     [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;
    eall_re=[eall_re; eps_re];
    tall_re=[tall_re ; t_re];
end
eall_re2=eall_re;
tall_re2=tall_re;

%~~
testnum=3
clear eall_re tall_re
eall_re=[];
tall_re=[];
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Nshift=size(REsamp.eps,3)
for whcase=1:Nshift    
    clear epsint_re x_resamp eps_re t_re    
    dz_samp=nanmean(diff(REsamp.z));
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
%    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.timeall(whcase,:) ) ;    
     [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;
    eall_re=[eall_re; eps_re];
    tall_re=[tall_re ; t_re];
end
eall_re3=eall_re;
tall_re3=tall_re;

%~~
testnum=4
clear eall_re tall_re
eall_re=[];
tall_re=[];
clear REsamp
% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Nshift=size(REsamp.eps,3)
for whcase=1:Nshift
    clear epsint_re x_resamp eps_re t_re    
    dz_samp=nanmean(diff(REsamp.z));
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
%    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.timeall(whcase,:) ) ;
     [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;
    eall_re=[eall_re; eps_re];
    tall_re=[tall_re ; t_re];
end
eall_re4=eall_re;
tall_re4=tall_re;

%~~

% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])

%idtr=isin(xx2.yday,[REsamp.timeall(1) REsamp.timeall(end)]);
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);
zvec=REsamp.z;

% real epsilon
epsint=nansum(xx2.eps)*10*1026;
[eps,t]=SimpleBoxCar(epsint(:,idtr),t_avg,xx2.yday(idtr));

%
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,4);

%~~
axes(ax(1))
ht=semilogy(t,eps,'k','linewidth',3)
hold on
% semilogy(tall_re3(1,:),eall_re3,'.','color','b')
% h3=semilogy(tall_re3(1,:),nanmean(eall_re3),'linewidth',2)

semilogy(tall_re3,eall_re3,'.','color',0.7*[1 1 1])
h3=semilogy(tall_re3,nanmean(eall_re3),'linewidth',2)

% plot +/- 1 std
semilogy(tall_re3(1,:),nanmean(eall_re3)-std(eall_re3),'--','color',0.3*[1 1 1])
semilogy(tall_re3(1,:),nanmean(eall_re3)+std(eall_re3),'--','color',0.3*[1 1 1])
ylim(yl)
xtloff
grid on
title(['Tchain ' num2str(whmoor) ' \rho\int \epsilon dz, ' num2str(t_avg) ' day avg. '])
ht=SubplotLetter(['w=0.15 m/s'],0.1,0.9)
set(ht,'fontsize',14)
ylabel('\rho \int \epsilon dz','fontsize',16)

%~~
axes(ax(2))
ht=semilogy(t,eps,'k','linewidth',3)
hold on
semilogy(tall_re1(1,:),eall_re1,'.','color',0.7*[1 1 1])
h1=semilogy(tall_re1(1,:),nanmean(eall_re1),'linewidth',2)
% plot +/- 1 std
semilogy(tall_re1(1,:),nanmean(eall_re1)-std(eall_re1),'--','color',0.3*[1 1 1])
semilogy(tall_re1(1,:),nanmean(eall_re1)+std(eall_re1),'--','color',0.3*[1 1 1])
ylim(yl)
xtloff
%ytloff
grid on
ht=SubplotLetter(['w=0.25 m/s'],0.3,0.9)
ylabel(' \rho \int \epsilon dz','fontsize',16)
set(ht,'fontsize',14)

%~~
axes(ax(3))
ht=semilogy(t,eps,'k','linewidth',3)
hold on
semilogy(tall_re2,eall_re2,'.','color',0.7*[1 1 1])
h2=semilogy(tall_re2(1,:),nanmean(eall_re2),'linewidth',2)
ylim(yl)
grid on
xtloff
ylabel('\rho \int \epsilon dz','fontsize',16)
ht=SubplotLetter(['w=0.5 m/s'],0.3,0.9)
set(ht,'fontsize',14)

%~~
axes(ax(4))
ht=semilogy(t,eps,'k','linewidth',3)
hold on
semilogy(tall_re4(1,:),eall_re4,'.','color',0.7*[1 1 1])
h4=semilogy(tall_re4(1,:),nanmean(eall_re4),'linewidth',2)
% plot +/- 1 std
semilogy(tall_re4(1,:),nanmean(eall_re4)-std(eall_re4),'--','color',0.3*[1 1 1])
semilogy(tall_re4(1,:),nanmean(eall_re4)+std(eall_re4),'--','color',0.3*[1 1 1])
ylim(yl)
%set(gca,'xtick',[170:5:185])
grid on
xlabel('Yearday','fontsize',16)
ylabel(' \rho \int \epsilon dz','fontsize',16)
ht=SubplotLetter(['w=0.75 m/s'],0.3,0.9)
set(ht,'fontsize',14)

shg

linkaxes(ax)

if saveplot==1    
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    ChkMkDir(figdir)
    figname=['Tchain' num2str(whmoor) '_Resamp_IntEpsvsTime_DiffSpeeds_2X2' ]
    print('-dpng',fullfile(figdir,figname))
end


%%