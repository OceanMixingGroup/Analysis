%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ResampleT2.m
%
% Resample 'T'2 mooring from T-tide to investigate effects of profiler
% sampling on overturns.
%
% Based on ResampleProfiler.m (for IWISE T-chains)
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
% 28 Feb. 2015  A. Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all


savedata=1

minOT=10 % minimum overturn size

useOrig=0
useGridded=1 ; dzgrid=10 ; dtgrid=2

addpath /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/

for  whtest=9%2:4
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
    elseif whtest==9
        w_samp=0.05; testnum=9 % m/s
    else
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        
    %~~  Load T-chain data to resample
    disp('loading data to resample')
    
    if useOrig==1
        source=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '.mat']
        load(source)
    elseif useGridded==1
        DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
        fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
        source=fullfile(DatDir,fname);
        load(source)
        
    else
    end
    
    % 'true' data to sample
    t_real=xx2.yday;
    z_real=xx2.z;
    data_real=xx2.T;
    
    % Sampling parameters
    clear SP
    SP.z_range=[xx2.z(1) xx2.z(end)];
    SP.dz_samp=10;
    SP.w_samp=w_samp;
    SP.t_start=18.5;     % start time
    SP.tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
    SP.time_range=2.5;
    
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
    REsamp=AddOverturnsToREsamp(REsamp,Params)
    
    
    if savedata==1
        DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
        ChkMkDir(DatDir)
        if useOrig==1
            fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
        elseif useGridded==1
            fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
        end
        save(fname,'REsamp')
    end
    
end


%%