%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_OT_LES.m
%
% Compute overturns with LES data at N2.
%
%
% * add temp/sgth to structure for later plotting and analysis
%
% 18 Feb. 2015
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

minOT=50;

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/

addpath /Users/Andy/Cruises_Research/Overturns_AP/

addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/

%load('N2_depth_time.mat')
load('N2_depth_time_Even.mat'); N2=N2even;clear N2even

[M,N]=size(N2.sgth);

d=nan*ones(size(N2.sgth));
eps=nan*ones(M,N);
Lot=nan*ones(M,N);
Lt=nan*ones(M,N);
Otnsq=nan*ones(M,N);

% overturns parameters
Params.lat=20.5;
Params.usetemp=0;
Params.minotsize=minOT;
Params.sigma=1e-5;
Params.runlmin=0;
Params.plotit=0;
P=sw_pres(N2.z,20.5);

hb=waitbar(0,'working')
for ind=1:N
    waitbar(ind/N,hb)
    [epsout,Lmin,Lotout,runlmax,Ltout,n2out,Otnsq_out,OT]=compute_overturns_discrete_AP_pden(N2.sgth(:,ind),P,Params);
    eps(:,ind)=epsout;
    Lot(:,ind)=Lotout;
    Lt(:,ind)=Ltout;
    Otnsq(:,ind)=Otnsq_out;
end
delete(hb)

%% plot epsilon and isopycnals

saveplot=1

thm=nanmean(N2.sgth,2);
dc=15;

close all
figure(1);clf
agutwocolumn(0.6)
wysiwyg
ezpc(N2.yday,N2.z,log10(eps));caxis([-9 -4.5]);cmap=jet;colormap([ 0.9*[1 1 1] ; cmap ])
%ezpc(N2.yday,N2.z,Lt);caxis([0 100]);colormap(jet)
hold on
contour(N2.yday,N2.z,N2.sgth,thm(1:dc:end),'k')
ylabel('Depth (m)')
cb=colorbar
cb.Label.String='log_{10} \epsilon'
xlabel('Yearday')
title('LES Model at N2')


if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['LES_N2_OTeps'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-png')
end

%% save results in structure 'OT'

OT=struct();
OT.sgth=N2.sgth;
OT.Lot=Lot;
OT.Lt=Lt;
OT.eps=eps;
OT.Otnsq=Otnsq;
OT.z=N2.z;
OT.yday=N2.yday
OT.Params=Params;
OT.MakeInfo=['Made ' datestr(now) ' w/ Compute_OT_LES.m in ' version]

fname=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_OT_minOT_' num2str(minOT) '.mat']
save(fname,'OT')

%%

ig=find(log10(eps)>-10);
figure(2);clf
%histogram(log10(eps(ig)));xlabel('log_{10}\epsilon');%xlim([-9 -5])
histogram(Lt(:));xlabel('L_T')
histogram(Lot(:));xlabel('L_{ot}')
%histogram(Otnsq(:));xlabel('L_{ot}')
%%