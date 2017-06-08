%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ComputeOverturnsTchainsReal.m
%
% Recompute overturns for JN's IWISE T-chain data, for use in resampling
% analyses. The T-chain data had an overturns field already, but something
% wasn't right so I recompute.
%
%
% 21 Nov. 2014 - AP - apickering@apl.washington.edu
% 26 Nov. - return overturn sizes also
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

whmoor=2

minOT=50;

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/

addpath /Users/Andy/Cruises_Research/Overturns_AP/

%~~  Load T-chain data to resample
load('all_moorings.mat')
load('all_gridded.mat')
c=[whmoor];
xx=grd{c};
clear grd

yday0=datenum2yday(xx.time);
if whmoor==3
    idt=find(yday0< 217);
elseif whmoor==4
    idt=find(yday0< 210);
elseif whmoor==1
    idt=find(yday0< 217);
elseif whmoor==2
    idt=find(yday0< 210);
else
    idt=1:length(xx.time);
end

xx2=xx;
xx2.time=xx.time(idt);
xx2.yday=datenum2yday(xx.time(idt));
xx2.T=xx.T(:,idt);
xx2.S=34.604-.045*(xx2.T-2.5);
xx2.eps=NaN*xx2.S;
xx2.Lot=NaN*xx2.S;
xx2.Lttot=NaN*xx2.S;
xx2.Lmin=NaN*xx2.S;
xx2.runlmax=NaN*xx2.S;
xx2.Otnsq_out=NaN*xx2.S;
xx2.d=NaN*xx2.S;
xx2.n2=NaN*xx2.S;

xx2.Lt_each=NaN*xx2.S;
xx2.Lot_each=NaN*xx2.S;
xx2.eps_each=NaN*xx2.S;
xx2.Otnsq_each=NaN*xx2.S;
%
hb=waitbar(0,'Wait a minute!')

% overturns parameters
lat=20.5;
usetemp=0;
minotsize=minOT;
sigma=1e-5;
runlmin=0;

Params.lat=lat;
Params.usetemp=usetemp;
Params.minotsize=minotsize;
Params.sigma=sigma;
Params.runlmin=runlmin;
Params.plotit=0;

for ind=1:length(xx2.time)
    waitbar(ind/length(xx2.time),hb)
    if mean(xx2.T(:,ind)<6)
        clear Epsout Lmin Lot runlmax Lttot n2
        % [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,minOT,1e-5,0);
        [Epsout,Lmin,Lot,runlmax,Lttot,ptmp,n2,Otnsq_out,OT]=compute_overturns_discrete_AP(xx2.z',xx2.T(:,ind),xx2.S(:,ind),Params);
        
        xx2.eps(:,ind)=Epsout;
        xx2.Lot(:,ind)=Lot;
        xx2.Lttot(:,ind)=Lttot;
        xx2.Lmin(:,ind)=Lmin;
        xx2.runlmax(:,ind)=runlmax;
        xx2.d(:,ind)=OT.d;
        xx2.n2(:,ind)=n2;
        xx2.Otnsq_out(:,ind)=Otnsq_out;
        
        clear NgoodOT
        NgoodOT=OT.Num_OT;%
        xx2.Lt_each(1:NgoodOT,ind)=OT.Lt_each(:);
        xx2.Lot_each(1:NgoodOT,ind)=OT.Lot_each(:);
        xx2.Otnsq_each(1:NgoodOT,ind)=OT.Otnsq_each(:);
        xx2.eps_each(1:NgoodOT,ind)=OT.eps_each(:);
        
    end
end

delete(hb)

% save data
xx2.OTparams=Params;
xx2.OTcode='compute_overturns_discrete_AP.m';
xx2.minOT=minOT;
xx2.DataSource='allgrided.mat';
xx2.MakeInfo=['Made ' datestr(now) ' w/ ComputeOverturnsTchainsReal.m']

DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) ]
ChkMkDir(DatDir)
fname=['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]

%
save(fullfile(DatDir,fname),'xx2')
%%