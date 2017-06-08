%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ResampleProfiler_V2.m
%
% Resample temperature from JN's IWISE (2011) T-chain moorings to simulate
% sampling by a MP or CTD.
%
% Newer version of ResampleProfiler.m. Now use more generalized function
% ResampleFieldGeneral.m and make REsamp structure directly instead of
% saving each case and combining after.
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
% Calls function ResampleFieldGeneral.m,AddOverturnsToREsamp.m
% 
%
%------------------------
% 17 Mar 2015 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

addpath /Users/Andy/Cruises_Research/SimProfiler/

savedata=1;
minOT=50; % minimum overturn size

for whmoor=3 % mooring #
    
    for  whtest=6
        clear testnum w_samp
        if whtest==1
            w_samp=0.25; testnum=whtest % m/s
        elseif whtest==2
            w_samp=0.5; testnum=whtest % m/s
        elseif whtest==3
            w_samp=0.15; testnum=whtest % m/s
        elseif whtest==4
            w_samp=0.75; testnum=whtest % m/s
        elseif whtest==5
            w_samp=1.0; testnum=whtest % m/s
        elseif whtest==6
            w_samp=1.5; testnum=whtest % m/s
        elseif whtest==7
            w_samp=2; testnum=whtest % m/s
        else
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        if whmoor==1
            TotalDays=50     % total # of days to sample
        elseif whmoor==2
            TotalDays=40     % total # of days to sample
        elseif whmoor==3
            %            TotalDays=50     % total # of days to sample
            TotalDays=20     % total # of days to sample
        elseif whmoor==4
            TotalDays=42
        else
        end
        
        %~~  Load T-chain data to resample
        disp('loading data to resample')
        %~~ use JN's T-chains during IWISE
        clear c xx t_real z_real data_real data_resamp mooring
        load('all_moorings.mat')
        load('all_gridded.mat')
        c=[whmoor];
        xx=grd{c};
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
        clear SP
        SP.z_range=[xx.z(1) xx.z(end)];
        SP.dz_samp=10;
        SP.w_samp=w_samp;
        SP.t_start=165;     % start time
        SP.tshift=(1/60/24) % shift each case by a few minutes to create ensemble of all phases
        SP.time_range=TotalDays;
        
        %
        REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)
        
        %% add overturns
        minOT=50;
        Params.lat=20.5;
        Params.usetemp=0;
        Params.minotsize=minOT;
        Params.sigma=1e-5;
        Params.runlmin=0;
        Params.plotit=0;
        REsamp=AddOverturnsToREsamp(REsamp,Params)
        
        %%
        
        %    save data
        if savedata==1
            DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Test' num2str(testnum)];
            ChkMkDir(DatDir)
            fname=fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
            REsamp.fname=fname;
            save(fname,'REsamp')
        end
        %
    end % whtest
    
end  % whmoor

%%