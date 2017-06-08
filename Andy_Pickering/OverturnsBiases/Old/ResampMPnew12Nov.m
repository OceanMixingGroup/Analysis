%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ResampMPnew12Nov.m
%
% Work on a better resampling to simulate profiler sampling of a data field
%
% 12 Nov 2014 A.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

savedata=1

minOT=50 % minimum overturn size

for whmoor=3%[3 4]%[1 2 3 4] % mooring #
    
    t_start=165;     % start time
    
    tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
    
    TotalDays=20     % total # of days to sample
    
    %~~ vertical velocity of profiler
    %   w_samp=0.25; testnum=1 % m/s
    % w_samp=0.5; testnum=2 % m/s
    %w_samp=0.15; testnum=3 % m/s
    %w_samp=0.75; testnum=4 % m/s
    %    w_samp=1.0; testnum=5 % m/s
    %        w_samp=0.4; testnum=6 % m/s
    w_samp=0.3; testnum=7 % m/s
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %
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
    z_real=xx.z;
    data_real=xx.T; % T-chain temperature
    data_resamp=nan*ones(1,length(t_real));
    %
    
    % depth-range to sample (m)
    z_range=[xx.z(1) xx.z(end)];
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    makefig=0
    
    hb=waitbar(0,'resampling data');
    
    % cycle through each case, creating an ensemble of phases
    for whcase=1:100
        waitbar(whcase/100,hb)
        close all
        clear t_range  t_grd t_samp z_samp t0 z0 data_resamp
        
        % start time
        t0=t_start + (whcase-1)*tshift;
        
        % starting depth of profiler
        z0=z_range(1);
        
        % define dz_samp
        dz_samp=10;
        dt_samp=dz_samp/w_samp;
        
        %~~~~~~ Everything below calculated from above chosen parameters
        
        % time to do 1 profile
        T_samp=range(z_range)/w_samp; % sec
        T_samp_mins=T_samp/60;        % min
        T_samp_day=T_samp/86400;      %days
        
        % How many profiles to do (each up or down ==1 profile)
        Nprof=round(TotalDays/T_samp_day);
        
        % depth-vector for sampling
        aZvec=z_range(1):dz_samp:z_range(end);
        
        %~~ # of points in one profile
        Nz=length(aZvec);
        
        dt_samp_day=dt_samp/86400;
        
        %~~ make time vector of sampling for all profiles
        t_samp=t0 : dt_samp_day : t0 + dt_samp_day*((Nz*Nprof)-1 );
        
        % make empty arrays
        zmat=nan*ones(Nz,Nprof);
        tmat=nan*ones(Nz,Nprof);
        data_resamp=nan*ones(Nz,Nprof);
        t_grid=nan*ones(1,Nprof);
        
        % resample along each profile segment
        for whp=1:Nprof
            clear tvev idt T Z datai
            tvec=t_samp( 1+ (whp-1)*Nz : Nz*whp);
            idt=isin(t_real,[tvec(1) tvec(end)]);
            
            clear zvec
            if iseven(whp)
                zvec=fliplr(aZvec);
            else
                zvec=aZvec;
            end
            
            datai= interp2(t_real(idt),z_real,data_real(:,idt),tvec,aZvec);
            data_resamp(:,whp)=datai(:);
            t_grid(whp)=nanmean(tvec);
            
            zmat(:,whp)=zvec;
            tmat(:,whp)=tvec;
            
        end
        %
        idtall=isin(t_real,[t_grid(1) t_grid(end)]);
        
        if makefig==1
            figure(1);clf
            
            ax1=subplot(211);
            ezpc(t_real(idtall),z_real,data_real(:,idtall))
            hold on
            plot(tmat,zmat,'w-')
            
            ax2=subplot(212);
            ezpc(t_grid,aZvec,data_resamp)
            xlabel('Yearday')
            ylabel('Depth')
            linkaxes([ax1 ax2])
        end
        
        % put results in structure 'x_resamp'
        x_resamp.minOT=minOT;
        x_resamp.tsampall=tmat;
        x_resamp.zall=zmat;
        x_resamp.time=t_grid;
        x_resamp.z=aZvec;
        x_resamp.T=data_resamp;
        x_resamp.S=34.604-.045*(x_resamp.T-2.5);
        x_resamp.eps=NaN*x_resamp.S;
        x_resamp.Lot=NaN*x_resamp.S;
        x_resamp.Lttot=NaN*x_resamp.S;
        %
        addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/
        
        % now compute overturns for these periods!
        
        for a=1:Nprof
            ind=a;
            if mean(x_resamp.T(:,ind)<6)
                clear Epsout Lmin Lot runlmax Lttot
                %        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.P(:,ind),xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
                [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(x_resamp.z',x_resamp.T(:,ind),x_resamp.S(:,ind),35.8,0,minOT,1e-5,0);
                x_resamp.eps(:,ind)=Epsout;
                x_resamp.Lot(:,ind)=Lot;
                x_resamp.Lttot(:,ind)=Lttot;
            end
        end
        
        x_resamp.w_samp=w_samp;
        x_resamp.MakeInfo=['Made ' datestr(now) ' w/ ResampMPnew12Nov.m in ' version]
        
        %    save data here so we can run a bunch of cases and plot after        
        if savedata==1
            fname=fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]);
            save(fname,'x_resamp')
        end
        
    end % whcase
    delete(hb)
    
end % whmoor
%%
