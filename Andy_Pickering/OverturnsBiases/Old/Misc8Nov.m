%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
% Code to resample a gridded depth-time field approximating the sampling of
% a moored profiler or CTD that continuously profiles up and down at a set
% speed.
%
% look at 1 turbulent event : 194:195
%
% 8 Nov. 2014 - AP - apickering@apl.washington.edu
%
% 10 Nov. - 25 days took 17 mins
%
% case 1: 165:170 ~ 1 spring
% case 1a: 165:190
% case1b:177:187 % ~ 1spring
% case1c: 190:200
%218:2  28 % ~ 1spring ; too many gaps in this period
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

savedata=0

t0=194;casename='c'
tstart=tic
tshift=(20/60/24)

%~~  Load T-chain data to resample
disp('loading data to resample')
%~~ use JN's T-chains during IWISE
load('all_moorings.mat')
load('all_gridded.mat')
c=[3];
xx=grd{c}
clear grd
t_real=datenum2yday(xx.time);
z_real=xx.z;
data_real=xx.T;
data_resamp=nan*ones(1,length(t_real));
%~~
whcase=1
for t_start=t0%:tshift:t0+10*tshift
    close all
    clear t_range z_range w_samp t_grd t_samp z_samp t0 z0 data_resamp
    
    t_range=[t_start t_start+2] % ~time range to sample (days)
    z_range=[580 2080] % depth-range to sample (m)
    
    %~~ vertical velocity of profiler
    w_samp=0.25 % m/s
    
    %~~ starting depth and time of profiler
    t0=t_range(1)
    z0=z_range(1);
    
    % define dt_sampe
%    dt_samp=8 % sampling dt (sec)
  %  dz_samp=w_samp*dt_samp % sampling dz (m);
  
  % define dz_samp
    dz_samp=2
    dt_samp=dz_samp/w_samp
    %~~~~~~ Everything below calculated from above chosen parameters
    
    % time to do 1 profile
    T_samp=range(z_range)/w_samp; % sec
    T_samp_mins=T_samp/60; % min
    T_samp_day=T_samp/86400; %days
    
    %~~ # of profiles to do (1 up or down=1 profile)
    % only works for even # profiles??
    Nprof=round(range(t_range)/T_samp_day);
    if ~iseven(Nprof);Nprof=Nprof+1;end
    
    % time of 'gridded MP' profiles==center of each profile
    % each angled (t-z) profile will get gridded into one 'vertical' profile at
    % the center time
    t_grid0=t0+T_samp_day/2
    t_grid=t_grid0  : T_samp_day : t_grid0+Nprof*T_samp_day    ;
    
    % begin and end times of each profile
    t_beg=nan*ones(Nprof,1);
    t_end=nan*ones(Nprof,1);
    t_beg=t0: T_samp_day : t0 + Nprof *T_samp_day;
    t_end=t0+T_samp_day :T_samp_day : t0 + Nprof*T_samp_day +T_samp_day;
    
    z_samp=[];
    nt_per_prof=T_samp/dt_samp;
    nt=nt_per_prof*Nprof;
    
    %~~ compute depth-time position of sampling. For 1st version, just do crude
    %loop. In future, could be simplified since z vector of each profile can
    %just be the same, flipped up or down...
    z_samp=[z0];
    t_samp=[t0];
    
    whdir=1;
    disp('Calculating sampling path')
    for wht=2:nt
        
        znew=z_samp(wht-1)+dz_samp;
        if znew>=z_range(2)
            whdir=-1;
        end
        
        if znew<=z_range(1)
            whdir=1;
        end
        
        z_samp=[z_samp z_samp(end)+whdir*dz_samp];
        t_samp=[t_samp t_samp(end)+dt_samp/86400];
        
    end
    
    %~~ Plot sampling path
    figure(1);clf
    plot(t_samp,z_samp,'o-')
    hold on
    axis ij
    ylabel('Depth (m) ')
    ylabel('Time (yearday)')
    
    % now try resampling a densly sampled (in time) data set
    % crudely  'resample' by finding real data points closest to each sampling point in
    % depth/time.  In future, maybe use interp2?
    disp('resampling data')
    hb=waitbar(0,'resampling data')
    for wht=1:length(z_samp)
        waitbar(wht/length(z_samp),hb)
        clear Iz It valz
        [valz,Iz]=nanmin(abs(z_samp(wht)-z_real));
        [valt,It]=nanmin(abs(t_samp(wht)-t_real));
        
        if valz<10 && valt<2/60/24
            data_resamp(wht)=data_real(Iz,It);
            %        figure(1)
            %      plot(t_real(It),z_real(Iz),'x')
        else
            disp('adfl')
            
        end
    end
    delete(hb)
    %
    data_resamp=data_resamp(:);
    data_gridded=nan*ones(nt_per_prof,Nprof);
    
    for whprof=1:Nprof-1
        clear idt
        idt=isin(t_samp,[t_beg(whprof) t_end(whprof)]);
        
        if nanmean(diff(z_samp(idt)))<0
            data_gridded(:,whprof)=flipud(data_resamp(idt));
        else
            data_gridded(:,whprof)=data_resamp(idt);
        end
    end
    
    %%
    cl=[1 8]
    yl=z_range
    figure(3);clf
    agutwocolumn(0.8)
    wysiwyg
    ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,2);
    
    axes(ax(1));
    ezpc(t_real,z_real,data_real)
    hold on
    plot(t_samp,z_samp,'w-','linewidth',1)
    ylim(yl)
    xlim([t_beg(1) t_end(end)])
    caxis(cl)
    colorbar
    
    axes(ax(2))
    %subplot(212)
    ezpc(t_grid(1:end-1),z_samp(idt),data_gridded)
    hold on
    plot(t_grid,z_range(1),'ko')
    ylim(yl)
    caxis(cl)
    colorbar
    xlabel('yearday')
    ylabel('Depth')
    
    linkaxes(ax)
    
    %    save2pdfAP('Tchain3_Resamp_164_Tcontour')
    %%
    
    % also get Get 'real' data profiles at times of center of MP profiles to
    % compare
    
    prof_real=nan*ones(length(z_real),length(t_grid));
    for whp=1:length(t_grid)
        [val,It]=nanmin(abs(t_grid(whp)-t_real));
        prof_real(:,whp)=data_real(:,It);
    end
    
    %% plot comparisons of 'gridded' profiles 1 at a time
    
    wht=1
    
    [val,It]=nanmin(abs(t_grid(wht)-t_real));
    figure(33);clf
    plot(data_real(:,[It]),z_real,'k','linewidth',2)
    hold on
    %plot(data_real(:,[It-30:2:It+30]),z_real,'color',0.5*[1 1 1]) % +/- 30 mins
    hold on
    
    plot(data_gridded(:,wht),z_samp(idt),'r','linewidth',2)
    axis ij
    
    %
    %% Now compute thorpe scales on resampled data
        
    %~~ T-chain dz is 10 m
    minOT=10; % min size of overturn
    
    %~~ First do resampled data
    
    x_resamp.time=t_grid;
    x_resamp.z=z_samp(idt);
    x_resamp.T=data_gridded;
    x_resamp.S=34.604-.045*(x_resamp.T-2.5);
    x_resamp.eps=NaN*x_resamp.S;
    x_resamp.Lot=NaN*x_resamp.S;
    x_resamp.Lttot=NaN*x_resamp.S;
    %
    addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/
    
    h = waitbar(0,'Please wait...');
    
    for a=1:length(x_resamp.time)-1
        ind=a;
        if mean(x_resamp.T(:,ind)<6)
            clear Epsout Lmin Lot runlmax Lttot
            %        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.P(:,ind),xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
            [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(x_resamp.z',x_resamp.T(:,ind),x_resamp.S(:,ind),35.8,0,minOT,1e-5,0);
            x_resamp.eps(:,ind)=Epsout;
            x_resamp.Lot(:,ind)=Lot;
            x_resamp.Lttot(:,ind)=Lttot;
        end
        %    if mod(a,1000)==0
        waitbar(a/length(x_resamp.time),h)
        %  end
    end
    
    delete(h)
    
    
    %% Plot results
    
    xx.yday=datenum2yday(xx.time);
    idtr=isin(xx.yday,[t_beg(1) t_end(end)]);
    %%
    
    %[eps_out,time_out]=SimpleBoxCar(nanmean(xx.eps),0.5,xx.yday);
    
    tavg=0.25
    
    dt_real=nanmean(diff(xx.yday))
    ntreal=round(tavg/dt_real)
    
    dtre=nanmean(diff(x_resamp.time));
    ntre=round(tavg/dtre)
    
    %~~
    figure(4);clf
    agutwocolumn(1)
    wysiwyg
    ax = MySubplot2(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,3);
    
    % plot timeseries of depth integrated eps
    axes(ax(1))
    semilogy(xx.yday(idtr),smooth(nanmean(xx.eps(:,idtr)),ntreal),'r')
    hold on
    hold on
    semilogy(x_resamp.time(1:end-1),smooth(nanmean(x_resamp.eps),ntre))
    
    cb=colorbar;killcolorbar(cb)
    grid on
    xtloff
    ylim([1e-8 1e-4])
    ht=SubplotLetter('6hr<\epsilon>')
    set(ht,'fontsize',15)
    
    legend('real','resamp','orientation','horziontal','location','best')

    % real epsilon, all profiles
    axes(ax(2));
    ezpc(xx.yday(idtr),xx.z,log10(xx.eps(:,idtr)))
    hold on
    ylabel('Depth ')
    xlim([t_grid(1) t_grid(end)])
    cb=    colorbar
    ylabel(cb,'log_{10}\epsilon')
    caxis([-9 -2])
    xtloff
    SubplotLetterMW('actual')
        
    % 'resampled' epsilon
    axes(ax(3));
    ezpc(x_resamp.time(1:end-1),x_resamp.z,log10(x_resamp.eps))
    colorbar
    caxis([-9 -2])
    xlabel('Yearday')
    cmap=jet;
    cmap=[0.7*[1 1 1];cmap;]
    colormap(cmap)
    SubplotLetterMW('resamp')
    
    linkaxes(ax,'x')
    
    %%
    %save2pdfAP('Tchain3_Resamp_164_171_epscontour')
    
    %% Plot time-average profiles
    
    figure(5);clf    
    semilogx(nanmean(x_resamp.eps,2),x_resamp.z,'linewidth',2)
    hold on
    semilogx(nanmean(xx.eps(:,idtr),2),xx.z,'m','linewidth',2)
    axis ij    
    grid on
    xlim([1e-8 1e-5])
    legend('resamp','real','location','best')
    ylabel('Depth')
    xlabel('<\epsilon>')
    title(['yday range= ' num2str(t_range(1)) ' - ' num2str(t_range(2)) ' '])
    %%
    %save2pdfAP('Tchain3_Resamp_164_171_epsProfile')
    
    if savedata==1
    fname=['Tchain3_resamp_case' num2str(whcase) casename]
    save(fname,'x_resamp')
    end
    whcase=whcase+1
end

telapsed = toc(tstart);
%%
    fname=['Tchain3_day194_196resamp']
    save(fname,'x_resamp')
    %%