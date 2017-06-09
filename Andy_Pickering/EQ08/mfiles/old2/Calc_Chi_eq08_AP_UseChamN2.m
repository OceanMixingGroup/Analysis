%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Calc_Chi_eq08_AP_UseChamN2.m
%
% Apply chipod method to EQ08 Chameleon profiles. In this version, use N2
% and dTdz computed in already-processed Chameleon data, instead of
% computing it myself from T and S. Want to check if differences are due to
% how I compute N2 and dT/dz.
%
%
% Modified from Calc_Chi_eq08_AP.m
%
%------------------
% 03/18/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

makeplots=0

Params.z_smooth=10;  % distance to smooth N^2 and dTdz over
Params.nfft=128;
Params.fmax=15;
Params.TPthresh=1e-6;
Params.resp_corr=1;  % correct TP spectra for freq response of thermistor
Params.gamma=0.2     % mixing efficiency

% choose directory to save processed casts in
datdirsave=fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/',...
    ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]);

% check if directory exists, make new if not
ChkMkDir(datdirsave)

% folder where processed cham files are saved
path_save='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/casts'

% make list of processed cham files
Flist=dir(fullfile(path_save,'*.mat'))
%%
% initialize a waitbar
hb=waitbar(0)

for icast=1:length(Flist)

    waitbar(icast/length(Flist),hb)
       
        
    try
        
        close all
        clear cal head
        % load the data for this cast
         load(fullfile(path_save,Flist(icast).name))
     
        clear cal
        cal=cal2;
        
        % average temp and sal in 1m bins like we normally do for CTD data
        clear zmin dz zmax tbin zbin sbin
        zmin=nanmin(cal.P);
        dz=1;
        zmax=nanmax(cal.P);
        minobs=2;
        [tbin zbin Nobs] = binprofile(cal.T ,cal.P, zmin, dz, zmax,minobs);
        [sbin zbin Nobs] = binprofile(cal.S,cal.P, zmin, dz, zmax,minobs);
        clear zmin dz zmax minobs
        
        %% compute dT/dz and N^2 for chi calculations
        clear ctd z_smooth
        ctd=struct();
        ctd.t1=tbin;
        ctd.s1=sbin;
        ctd.p=zbin;
        ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        ctd=Compute_N2_dTdz_forChi(ctd,Params.z_smooth);
        
        % are pressure loops already removed? or maybe there are none/few?
        
        % make windows for chi calcs
        clear TP
        TP=cal.TP;
        [todo_inds,Nwindows]=MakeCtdChiWindows(TP,Params.nfft);
        %    [todo_inds,Nwindows]=MakeCtdChiWindows_simple(TP,Params.nfft);
        
        %~ make 'avg' structure for the processed data
        clear avg
        avg=struct();
        avg.Params=Params;
        tfields={'datenum','P','N2','dTdz','fspd','T','S','P','theta','sigma',...
            'chi1','eps1','KT1','TP1var','kstart','kstop'};
        for n=1:length(tfields)
            avg.(tfields{n})=NaN*ones(Nwindows,1);
        end
        
        % get TP samplerate (Hz)
        avg.samplerate=head.slow_samp_rate * head.irep.TP;
        
        % get average time, pressure, and fallspeed in each window
        for iwind=1:Nwindows
            clear inds
            inds=todo_inds(iwind,1) : todo_inds(iwind,2);
            avg.P(iwind)=nanmean(cal.P(inds));
            avg.fspd(iwind)=nanmean(cal.fspd(inds))/100;
        end
        
        
        %~ 02/17/16 - AP - try saving spectra to plot after
        %     tpspec=nan*ones(Nwindows,Params.nfft /2);
        %     kspec=nan*ones(Nwindows,Params.nfft /2);
        %     kkspec=nan*ones(Nwindows,56);
        %~
        
        
        if makeplots==1
            %~~ plot histogram of avg.P to see how many good windows we have in
            %each 10m bin
            if nanmax(avg.P>10)
                figure
                hi=histogram(avg.P,0:10:nanmax(avg.P));
                hi.Orientation='Horizontal';axis ij;
                ylabel('P [db]')
                xlabel('# good data windows')
                % title([whSN ' cast ' cast_suffix ' - ' C.castdir 'cast'],'interpreter','none')
                % print('-dpng',fullfile(chi_fig_path_specific,[whSN '_' cast_suffix '_Fig' num2str(whfig) '_' C.castdir 'cast_chi_' whsens '_avgPhist']))
                % whfig=whfig+1;
            end
        end
        
        
        % get N2, dTdz for each window
        good_inds=find(~isnan(ctd.p));
        % interpolate ctd data to same pressures as chipod
        avg.N2  =interp1(ctd.p(good_inds),ctd.N2(good_inds),avg.P);
        avg.dTdz=interp1(ctd.p(good_inds),ctd.dTdz(good_inds),avg.P);
        avg.T   =interp1(ctd.p(good_inds),ctd.t1(good_inds),avg.P);
        avg.S   =interp1(ctd.p(good_inds),ctd.s1(good_inds),avg.P);
        
        % note sw_visc not included in newer versions of sw?
        avg.nu=sw_visc_ctdchi(avg.S,avg.T,avg.P);
        avg.tdif=sw_tdif_ctdchi(avg.S,avg.T,avg.P);
        
        
        % loop through each window and do the chi computation
        for iwind=1:Nwindows
            
            clear inds
            % get inds for this window
            inds=todo_inds(iwind,1) : todo_inds(iwind,2);
            
            % integrate dT/dt spectrum
            clear tp_power freq
            [tp_power,freq]=fast_psd(TP(inds),Params.nfft,avg.samplerate);
            avg.TP1var(iwind)=sum(tp_power)*nanmean(diff(freq));
            
            if avg.TP1var(iwind)>Params.TPthresh
                
                % apply filter correction for sensor response?
                
                %figure(10);clf
                %loglog(freq,tp_power)
                
                
                if Params.resp_corr==1
                    trans_fcn=0;
                    trans_fcn1=0;
                    thermistor_filter_order=2;
                    thermistor_cutoff_frequency=32;
                    analog_filter_order=4;
                    analog_filter_freq=50;
                    tp_power=invert_filt(freq,invert_filt(freq,tp_power,thermistor_filter_order, ...
                        thermistor_cutoff_frequency),analog_filter_order,analog_filter_freq);
                end
                
                %             hold on
                %             loglog(freq,tp_power)
                %             legend('raw','corr')
                %             xlim([1e0 1e2])
                %             freqline(10)
                %             freqline(20)
                %             grid on
                %             pause
                
                % compute chi using iterative procedure
                [chi1,epsil1,k,spec,kk,speck,stats]=get_chipod_chi(freq,tp_power,abs(avg.fspd(iwind)),avg.nu(iwind),...
                    avg.tdif(iwind),avg.dTdz(iwind),'nsqr',avg.N2(iwind),'fmax',Params.fmax,'gamma',Params.gamma);%,'doplots',1);
                %            pause
                %            'doplots',1 for plots
                avg.chi1(iwind)=chi1(1);
                avg.eps1(iwind)=epsil1(1);
                avg.KT1(iwind)=0.5*chi1(1)/avg.dTdz(iwind)^2;
                avg.kstart(iwind)=stats.k_start;
                avg.kstop(iwind)=stats.k_stop;
                
                %             % 02/17/16 - AP - save spectra
                %             tpspec(iwind,:)=tp_power;
                %             kspec(iwind,:)=spec;
                %             kkspec(iwind,:)=speck;
                %             fspec=freq;
                %             ks=k;
                %
                %             if ~isnan(kk)
                %             kks=kk;
                %             end
                
            end % if T1Pvar>threshold
            
        end % windows
        
        
        if makeplots==1
            %~ plot summary figure
            ax=CTD_chipod_profile_summary(avg,cal,TP);
            axes(ax(1))
            title(Flist(icast).name,'interpreter','none')
            %~~~
        end
        
        %~
        %     avg.tpspec=tpspec;
        %     avg.kspec=kspec;
        %     avg.kkspec=kkspec;
        %     avg.ks=ks;
        %     avg.kks=kks;
        %     avg.fspec=fspec;
        
        %~
        
        avg.lat=ctd.lat;
        avg.z=sw_dpth(avg.P,avg.lat);
        avg.MakeInfo=['Made ' datestr(now) ' w/ Calc_Chi_eq08_AP.m'];
        
        % save results
        save(fullfile(datdirsave,[Flist(icast).name(1:end-4) 'avg.mat']),'avg')
        
    end % try
    
    
    
end % icast

delete(hb)


%%