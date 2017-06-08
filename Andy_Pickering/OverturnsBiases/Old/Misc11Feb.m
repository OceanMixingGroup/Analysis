%%
%
% Misc11Feb.m
%
%
%% figure out how overturns code works in detail 

clear ; close all

whmoor=4

minOT=50;

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/
    
%~~  Load T-chain data to resample
load('all_moorings.mat')
load('all_gridded.mat')
c=[whmoor];
xx=grd{c};
clear grd mooring config

%%
%t_real=datenum2yday(xx.time);
xx.yday=datenum2yday(xx.time);
%idtr=isin(xx.yday,[x_resamp.time(1) x_resamp.time(end)]);

xx2=xx;
xx2.time=xx.time;
xx2.yday=xx.yday;
xx2.T=xx.T;
xx2.S=34.604-.045*(xx2.T-2.5);
xx2.eps=NaN*xx2.S;
xx2.Lot=NaN*xx2.S;
xx2.Lttot=NaN*xx2.S;
xx2.Lmin=NaN*xx2.S;
xx2.runlmax=NaN*xx2.S;

ind=15100
mean(xx2.T(:,ind)<6)
clear Epsout Lmin Lot runlmax Lttot
[Epsout,Lmin,Lot,runlmax,Lttot,ptmp]=compute_overturns_discrete_APcomment(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,minOT,1e-5,0);
%[Epsout,Lmin,Lot,runlmax,Lttot,ptmp]=compute_overturns_discrete(p,t,s,lat,usetemp,minotsize,sigma,runlmin);
% p=xx2.z';
% t=xx2.T(:,ind);
% s=xx2.S(:,ind);
% 35.8,0,minOT,1e-5,0);

%        xx2.eps(:,ind)=Epsout;

figure(1);clf
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.06, 4,1);

axes(ax(1))
plot(xx2.T(:,ind),xx2.z,'o-')
hold on
plot(ptmp,xx2.z)
axis ij
ylim([1200 2200])
title(['ind=' num2str(ind) ' ,yday=' num2str(xx2.yday(ind))])
grid on
ylabel('Depth (m) ')
xlabel('T')

axes(ax(2))
plot(Lot,xx2.z,'o-')
axis ij
ylim([1200 2200])
%title(['ind=' num2str(ind) ' ,yday=' num2str(xx2.yday(ind))])
grid on
xlim([0 max(Lot)+50])
xlabel('L_T')
ytloff


axes(ax(3))
plot(Lttot,xx2.z,'o-')
axis ij
ylim([1200 2200])
%title(['ind=' num2str(ind) ' ,yday=' num2str(xx2.yday(ind))])
grid on
xlim([0 max(Lot)+50])
xlabel('L')
ytloff

axes(ax(4))
plot(log10(Epsout),xx2.z,'o-')
axis ij
ylim([1200 2200])
%title(['ind=' num2str(ind) ' ,yday=' num2str(xx2.yday(ind))])
grid on
xlim([-11 -4])
xlabel('log_{10} \epsilon')
ytloff


linkaxes(ax,'y')
%%