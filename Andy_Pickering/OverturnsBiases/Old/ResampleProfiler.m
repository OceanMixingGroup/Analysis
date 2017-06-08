%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ResampleProfiler.m
%
% ** NOTE see new version ResampleProfiler_V2.m **
%
% Resample temperature from JN's IWISE T-chain moorings to simulate
% sampling by a MP or CTD.
%
% Based on previous version: ResampMPnew12Nov.m . Modified to use
% MakeProfPath.m function to create sampling path instead of doing in loop.
%
%
% 6 Feb 2015 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

savedata=1

minOT=50 % minimum overturn size

for whmoor=1 % mooring #
    
    for  whtest=1:4
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
        elseif whtest==8
            w_samp=0.1; testnum=8 % m/s
        else
        end
        
        %    w_samp=0.4; testnum=6 % m/s
        %    w_samp=0.3; testnum=7 % m/s
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        t_start=165;     % start time
        
        tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
        
        if whmoor==1
            TotalDays=50     % total # of days to sample
        elseif whmoor==2
            TotalDays=40     % total # of days to sample
        elseif whmoor==3
            TotalDays=50     % total # of days to sample
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
        
        % depth-range to sample (m)
        clear z_range dz_samp dztot Tprof Nshift
        z_range=[xx.z(1) xx.z(end)];
        dz_samp=10;
        
        dztot=range(z_range);
        Tprof=dztot/w_samp/86400; % time for a profile (1 up or down) in days
        Nshift=round(2*Tprof/tshift)% number of shifts to go through all phases
        
        makefig=0;
        hb=waitbar(0,'resampling data');
        
        % cycle through each case, creating an ensemble of phases
        for whcase=1:Nshift
            waitbar(whcase/Nshift,hb)
            %close all
            
            clear t0 time_range tmat zmat Nprof Nz
            t0=t_start + (whcase-1)*tshift;
            time_range=TotalDays;
            plotit=0;
            [tmat,zmat]=MakeProfPath(t0,z_range,dz_samp,time_range,w_samp,plotit);
            
            Nprof=size(zmat,2);
            Nz=size(zmat,1);
            
            % make empty arrays
            clear data_resamp t_grid
            data_resamp=nan*ones(Nz,Nprof);
            t_grid=nan*ones(1,Nprof);
            
            % resample along each profile segment
            for whp=1:Nprof
                clear tvec idt datai
                tvec=tmat(:,whp);
                idt=isin(t_real,[tvec(1) tvec(end)]);
                %                datai= interp2(t_real(idt),z_real,data_real(:,idt),tvec,zmat(:,1));
                datai= interp2(t_real(idt),z_real,data_real(:,idt),tvec,zmat(:,whp));
                if nanmean(diff(zmat(:,whp)))<0
                    data_resamp(:,whp)=flipud(datai(:));
                else
                    data_resamp(:,whp)=datai(:);
                end
                %                data_resamp(:,whp)=datai(:);
                t_grid(whp)=nanmean(tvec);
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
                ezpc(t_grid,zmat(:,1),data_resamp)
                xlabel('Yearday')
                ylabel('Depth')
                linkaxes([ax1 ax2])
            end
            
            % put results in structure 'x_resamp'
            clear x_resamp
            x_resamp.minOT=minOT;
            x_resamp.tsampall=tmat;
            x_resamp.zall=zmat;
            x_resamp.time=t_grid;
            x_resamp.z=zmat(:,1);
            x_resamp.T=data_resamp;
            x_resamp.S=34.604-.045*(x_resamp.T-2.5);
            x_resamp.eps=NaN*x_resamp.S;
            x_resamp.Lot=NaN*x_resamp.S;
            x_resamp.Lttot=NaN*x_resamp.S;
            x_resamp.d=NaN*x_resamp.S;
            x_resamp.Otnsq_out=NaN*x_resamp.S;
            x_resamp.n2=NaN*x_resamp.S;
            
            x_resamp.Lt_each=NaN*x_resamp.S;
            x_resamp.Lot_each=NaN*x_resamp.S;
            x_resamp.Otnsq_each=NaN*x_resamp.S;
            x_resamp.eps_each=NaN*x_resamp.S;
            
            x_resamp.start=NaN*x_resamp.S;
            
            addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/
            
            %~ Now compute overturns for these periods!
            
            % overturns parameters
            Params.lat=20.5;
            Params.usetemp=0;
            Params.minotsize=minOT;
            Params.sigma=1e-5;
            Params.runlmin=0;
            Params.plotit=0;
            
            addpath /Users/Andy/Cruises_Research/Overturns_AP/
            
            for a=1:Nprof
                ind=a;
                if mean(x_resamp.T(:,ind)<6)
                    clear Epsout Lmin Lot runlmax Lttot OT Otnsq_out
                    %[Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.P(:,ind),xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
                    % [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(x_resamp.z(:),x_resamp.T(:,ind),x_resamp.S(:,ind),35.8,0,minOT,1e-5,0);
                    [Epsout,Lmin,Lot,runlmax,Lttot,ptmp,n2,Otnsq_out,OT]=compute_overturns_discrete_AP(x_resamp.z(:),x_resamp.T(:,ind),x_resamp.S(:,ind),Params);
                    x_resamp.eps(:,ind)=Epsout;
                    x_resamp.Lot(:,ind)=Lot;
                    x_resamp.Lttot(:,ind)=Lttot;
                    x_resamp.d(:,ind)=OT.d;
                    x_resamp.Otnsq_out(:,ind)=Otnsq_out;
                    x_resamp.n2(:,ind)=n2;
                    
                    clear NgoodOT
                    NgoodOT=OT.Num_OT;%
                    
                    x_resamp.Lt_each(1:NgoodOT,ind)=OT.Lt_each(:);
                    x_resamp.Lot_each(1:NgoodOT,ind)=OT.Lot_each(:);
                    x_resamp.Otnsq_each(1:NgoodOT,ind)=OT.Otnsq_each(:);
                    x_resamp.eps_each(1:NgoodOT,ind)=OT.eps_each(:);
                    
                end
            end % Nprof
            
            x_resamp.Nprof=Nprof;
            x_resamp.w_samp=w_samp;
            x_resamp.MakeInfo=['Made ' datestr(now) ' w/ ResampleProfiler.m in ' version];
            
            %    save data here so we can run a bunch of cases and plot after
            if savedata==1
                DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Test' num2str(testnum)];
                if whcase==1
                    ChkMkDir(DatDir)
                end
                fname=fullfile(DatDir,['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]);
                x_resamp.fname=fname;
                save(fname,'x_resamp')
            end
            
        end % whcase
        delete(hb)
        
    end % whtest
    
end  % whmoor

%%
