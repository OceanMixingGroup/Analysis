%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Resample_LES_OT.m
%
% Resample LES data at N2 to simulate profiling instrument at different
% speeds and assess effects on overturns.
%
% Based on ResampleProfiler.m, written for IWISE T-chain sampling analysis.
%
% OUTPUT
% Results are output in a structure 'REsamp' with fields:
%
% w_samp      : Vertical sampling speed of 'profiler'
% Nshift      : Number of sampling ensembles (each shifted slightly in time)
% data_resamp : The resampled data [Depth X time X Nshift]
% data_real   : the 'true' data to be resampled [Depth X time]
% tsamp       : Time of resampling [Depth X time X Nshift]
% zsamp       : Depths of resampling [Depth X time X Nshift]
% tgrid       : Mid-time of each sample profile (time that would be assigned to
% profiler data)
% zreal       : Depth vector of data_real
% treal       : Time-vector of data_real
% SP          : Sampling parameters for ResampleFieldGeneral.m

%
% 19 Feb. 2015 - A. Pickering - apickering@coas.oregonstate.edu
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES

addpath /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/

savedata=1

minOT=50 % minimum overturn size
makefig=0;

for  whtest=1:4 % various sampling speeds
    
    close all
    
    clear testnum w_samp
    if whtest==1
        w_samp=0.25; testnum=1 % m/s
    elseif whtest==2
        w_samp=0.5; testnum=2 % m/s
    elseif whtest==3
        w_samp=0.15; testnum=3 % m/s
    elseif whtest==4
        w_samp=0.75; testnum=4 % m/s
    elseif whtest==5
        w_samp=1.0; testnum=5 % m/s
    else
    end
    
    %~~  Load T-chain data to resample
    disp('loading data to resample')
    load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/N2_depth_time_Even.mat')
    N2=N2even;clear N2even;
    TotalDays=range(N2.yday); % total # of days to sample
    
    data_real=N2.sgth;
    t_real=N2.yday;
    z_real=N2.z;
    
    % Sampling parameters
    clear SP
    SP.z_range=[N2.z(1) N2.z(end)];
    SP.dz_samp=10;
    SP.w_samp=w_samp;
    SP.t_start=238.9602;     % start time
    SP.tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
    SP.time_range=TotalDays;
    %
    REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)
    
    
    %% add overturns
    minOT=minOT;
    Params.lat=20.5;
    Params.usetemp=0;
    Params.minotsize=minOT;
    Params.sigma=1e-5;
    Params.runlmin=0;
    Params.plotit=0;
    
    %    REsamp=AddOverturnsToREsamp(REsamp,Params)
    REsamp=AddOverturnsToREsamp_pden(REsamp,Params);
    %%
    %    save data
    if savedata==1
        DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/Data/Test' num2str(testnum)]
        ChkMkDir(DatDir)
        fname=fullfile(DatDir,['LES_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
        save(fname,'REsamp')
    end
    
    
end % whtest


%%
