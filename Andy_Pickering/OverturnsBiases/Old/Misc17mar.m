
clear ; close all ; clc

% Load data to resample

whmoor=3
%~~  Load T-chain data to resample
disp('loading data to resample')
%~~ use JN's T-chains during IWISE
clear c xx t_real z_real data_real data_resamp mooring
load('all_moorings.mat')
load('all_gridded.mat')
c=[whmoor];
xx=grd{c}
clear grd mooring config
t_real=datenum2yday(xx.time);

if whmoor==3
    idt=find(t_real< 217);
    t_real=t_real(idt);
    data_real=xx.T(:,idt);
    clear idt
elseif whmoor==4
    idt=find(t_real< 210);
    t_real=t_real(idt);
    data_real=xx.T(:,idt);
    clear idt
else
    data_real=xx.T; % T-chain temperature
end

z_real=xx.z;

clear ; close all ; clc

% Load data to resample

whmoor=3
%~~  Load T-chain data to resample
disp('loading data to resample')
%~~ use JN's T-chains during IWISE
clear c xx t_real z_real data_real data_resamp mooring
load('all_moorings.mat')
load('all_gridded.mat')
c=[whmoor];
xx=grd{c}
clear grd mooring config
t_real=datenum2yday(xx.time);

if whmoor==3
    idt=find(t_real< 217);
    t_real=t_real(idt);
    data_real=xx.T(:,idt);
    clear idt
elseif whmoor==4
    idt=find(t_real< 210);
    t_real=t_real(idt);
    data_real=xx.T(:,idt);
    clear idt
else
    data_real=xx.T; % T-chain temperature
end

z_real=xx.z;

% Sampling parameters
SP.w_samp=0.5
SP.tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
SP.t_start=165
SP.time_range=5;
SP.z_range=[xx.z(1) xx.z(end)];
SP.dz_samp=10;


%% resample data

REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)
%REsamp=ResampleFieldGeneral(data_real,t_real,z_real)

%% add overturns
REsamp=AddOverturnsToREsamp(REsamp)
%%
%
%
%%

clear ; close all


load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3/Test1/Tchain3_Test1_minOT_50_AllCases.mat')
RE1=REsamp;clear REsamp
load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3/Test1/Tchain3_Test1_minOT_50_AllCases_V2.mat')
RE2=REsamp;clear REsamp

%%
whc=5
figure(1);clf
plot(RE1.tsamp(:,:,whc),RE1.zsamp(:,:,whc),'k')
hold on
plot(RE2.tsamp(:,:,whc),RE2.zsamp(:,:,whc),'m--')
%%

dif=RE1.t(:,:,5)-RE2.data_resamp(:,:,5)
dif=RE1.eps(:,:,5)-RE2.eps(:,:,5)
%%