%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ComputeOverturnsT2_gridded.m
%
% Same as ComputeOverturnsT2.m, except use gridded data made in
% Make_T2_grid.m
%
% Compute overturns at T2. Compare results for raw, and interpolated to
% even depth grid.
%
% 2 Mar. 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

savedata=1

minOT=20;

dzgrid=20 % depth spacing in meters
dtgrid=2  % time spacing in minutes

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/

addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/

addpath /Users/Andy/Cruises_Research/Overturns_AP/

%~~  Load T-chain data to resample

Bdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/'
source=fullfile(Bdir,['T2_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min'])

%source='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_Gridded'
load(source)


xx2=struct();
xx2.z=T2.z;
%xx2.time=dat.time;
xx2.yday=T2.yday;%datenum2yday(xx2.time);
xx2.T=T2.t;
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

for ind=1:length(xx2.yday)
    waitbar(ind/length(xx2.yday),hb)
    %    if mean(xx2.T(:,ind)<6)
    numg=find(~isnan(xx2.T(:,ind)));
    if length(numg)>2
        clear Epsout Lmin Lot runlmax Lttot n2
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
xx2.DataSource=source;
xx2.MakeInfo=['Made ' datestr(now) ' w/ ComputeOverturnsT2_gridded.m']

DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
ChkMkDir(DatDir)
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']

%
if savedata==1
save(fullfile(DatDir,fname),'xx2')
end
%%