%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ComputeChi_Chameleon_tiwe.m
%
% Apply chi-pod processing methods to chameleon FP07 data from tiwe profiles.
%
% Chameleon mat files are made w/ Run_tiwe_AP.m, which must be run
% first.
%
% After this, combine all profiles w/ Combine_XXXX_chi.m
%
%
%------------------------
% 04/7/17 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

makeplots = 0
savespec  = 1  % option to save wavenumber spectra

% Params for chipod calculations
Params.z_smooth = 1;     % distance to smooth N^2 and dTdz over
Params.nfft     = 128;   %
Params.fmax     = 10;     %
Params.TPthresh = 1e-6;  %
Params.resp_corr= 0;     % correct TP spectra for freq response of thermistor
Params.fc       = 99     % cutoff frequency for response correction
Params.gamma    = 0.2    % mixing efficiency

% specify method of computing N^2 and dT/dz
whN2dTdz  = 'regular'
%whN2dTdz = 'regular2'
%whN2dTdz = 'line'
%whN2dTdz = 'raw_line'

% Add all the paths we need from mixing software
mixpath = '/Users/Andy/Cruises_Research/mixingsoftware/' ;
addpath(fullfile(mixpath,'seawater'))
addpath(fullfile(mixpath,'general'))
addpath(fullfile(mixpath,'marlcham'))
addpath(fullfile(mixpath,'CTD_Chipod','mfiles')) ;
addpath(fullfile(mixpath,'chipod','compute_chi')); % get_chipod_chi.md

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/code/

tiwe_patches_paths

%%
if Params.resp_corr==0
    Params.fc=99;
end

% Make directory to save processed casts in (name based on Params)
if strcmp(whN2dTdz, 'regular')
    datdirsave=fullfile(path_chipod_bin,...
        ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft)]);
else
    datdirsave=fullfile(path_chipod_bin,...
        ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft) '_' whN2dTdz]);
end

disp(['Data will be saved to ' datdirsave])

% check if directory exists, make new if not
ChkMkDir(datdirsave)

% initialize a waitbar
hb=waitbar(0);

%%
tstart=tic;
% loop through each cast
for icast=3870:3918
    
    % update waitbar
    waitbar(icast/3918,hb,['time elapsed=' num2str(round(toc(tstart)/60)) 'mins'])
    
    try
        
        close all
        clear cal head
        clear tpspec kspec kkspec fspec kks ks
        
        % Load the data for this cast
        load(fullfile(save_dir_cal,['tw91_' sprintf('%04d',icast) '_cal.mat']))
        clear cal
        cal     = cal2;
        
        % Average temp and sal in 1m bins like we normally do for CTD data
        clear zmin dz zmax tbin zbin sbin
        zmin=nanmin(cal.P);
        dz=1;
        zmax=nanmax(cal.P);
        minobs=2;
        [tbin zbin Nobs] = binprofile(cal.T1 ,cal.P, zmin, dz, zmax, minobs);
        [sbin zbin Nobs] = binprofile(cal.SAL,cal.P, zmin, dz, zmax, minobs);
        clear zmin dz zmax minobs
        
        % Put in 'ctd' structure
        clear ctd z_smooth
        ctd = struct();
        ctd.t1 = tbin;
        ctd.s1 = sbin;
        ctd.p  = zbin;
        
        % add in lat and lon (from sum_eq14.m)
%         clear idot lat1 lat2
%         idot = strfind(head.lon.start,'.');
%         lon1 = str2num(head.lon.start(1:idot-3));
%         lon2 = str2num(head.lon.start(idot-2:end))/60;
%         ctd.lon = nanmean([lon1 lon2]);
        
%         clear idot lat1 lat2
%         idot = strfind(head.lat.start,'.');
%         lat1 = str2num(head.lat.start(1:idot-3));
%         lat2 = str2num(head.lat.start(idot-2:end))/60;
%         ctd.lat = nanmean([lat1 lat2]);
ctd.lat=0.3;
        %    ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        
        % compute N^2 and dT/dz
        if strcmp(whN2dTdz, 'regular')
            ctd = Compute_N2_dTdz_forChi(ctd,Params.z_smooth);
        elseif strcmp(whN2dTdz, 'regular2')
            ctd = Compute_N2_dTdz_forChi_2(ctd);
        elseif strcmp(whN2dTdz, 'line')
            ctd = Compute_N2_dTdz_forChi_line(ctd,1);
        elseif strcmp(whN2dTdz, 'raw_line')
            % use raw (24hz) CTD data to comptue N^2 and Tz
            ctd = Compute_N2_dTdz_forChi_raw_line(ctd,cal,0.5);
        end
        
        % are pressure loops already removed? or maybe there are none/few?
        
        % Make windows for chi calcs
        clear TP
        TP = cal.TP;
        [todo_inds,Nwindows] = MakeCtdChiWindows(TP,Params.nfft);
        
        %~ make 'avg' structure for the processed data
        clear avg
        avg=struct();
        avg.Params=Params;
        tfields={'datenum','P','N2','dTdz','fspd','T','S','P','theta','sigma',...
            'chi1','eps1','KT1','TP1var','kstart','kstop'};
        for n=1:length(tfields)
            avg.(tfields{n})=NaN*ones(Nwindows,1);
        end
        
        avg.FspdStd = avg.fspd;
        
        % Get TP samplerate (Hz)
        avg.samplerate = head.slow_samp_rate * head.irep.T1P;
        
        % Get average time, pressure, and fallspeed in each window
        for iwind=1:Nwindows
            clear inds
            inds = todo_inds(iwind,1) : todo_inds(iwind,2);
            avg.P(iwind)    = nanmean(cal.P(inds));
            avg.fspd(iwind) = nanmean(cal.FALLSPD(inds))/100;
            avg.FspdStd(iwind) = nanstd(cal.FALLSPD(inds)/100);
        end
        
        %~ save spectra to plot after
        % make empty arrays for spectra
        if savespec==1
            % observed wavenumber
            ks     = nan*ones(Nwindows,Params.nfft /2);
            % observed wave# spectra
            kspec  = nan*ones(Nwindows,Params.nfft /2);
            % fit wavenumber
            kks    = nan*ones(Nwindows,56);
            % fit wave# spectra
            kkspec = nan*ones(Nwindows,56);
            % observed frequency spectra
            tpspec = nan*ones(Nwindows,Params.nfft /2);
            %~
        end
        
        
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
        
        
        %
        good_inds = find(~isnan(ctd.p));
        % interpolate ctd data to same pressures as chipod
        avg.N2  = interp1(ctd.p(good_inds), ctd.N2(good_inds)  ,avg.P );
        avg.dTdz= interp1(ctd.p(good_inds), ctd.dTdz(good_inds),avg.P );
        avg.T   = interp1(ctd.p(good_inds), ctd.t1(good_inds)  ,avg.P );
        avg.S   = interp1(ctd.p(good_inds), ctd.s1(good_inds)  ,avg.P );
        
        % note sw_visc not included in newer versions of sw?
        avg.nu   = sw_visc_ctdchi(avg.S,avg.T,avg.P);
        avg.tdif = sw_tdif_ctdchi(avg.S,avg.T,avg.P);
        
        avg.n_iter = nan*ones(size(avg.P));
        
        % loop through each window and do the chi computation
        for iwind=1:Nwindows
            
            clear inds
            % get inds for this window
            inds = todo_inds(iwind,1) : todo_inds(iwind,2);
            
            % integrate dT/dt spectrum
            clear tp_power freq
            [tp_power,freq] = fast_psd(TP(inds),Params.nfft,avg.samplerate);
            avg.TP1var(iwind) = sum(tp_power)*nanmean(diff(freq));
            
            if avg.TP1var(iwind)>Params.TPthresh
                
                % apply filter correction for sensor response?
                
                if Params.resp_corr==1
                    thermistor_filter_order=2; % 2=double pole?
                    thermistor_cutoff_frequency=Params.fc;%32;
                    analog_filter_order=4;
                    analog_filter_freq=50;
                    tp_power=invert_filt(freq,invert_filt(freq,tp_power,thermistor_filter_order, ...
                        thermistor_cutoff_frequency),analog_filter_order,analog_filter_freq);
                end
                
                % Compute chi using iterative procedure
                [chi1,epsil1,k,spec,kk,speck,stats] = get_chipod_chi(freq,tp_power,abs(avg.fspd(iwind)),avg.nu(iwind),...
                    avg.tdif(iwind),avg.dTdz(iwind),'nsqr',avg.N2(iwind),'fmax',Params.fmax,'gamma',Params.gamma);%,'doplots',1);
                %            pause
                %            'doplots',1 for plots
                avg.chi1(iwind)   = chi1(1);
                avg.eps1(iwind)   = epsil1(1);
                avg.KT1(iwind)    = 0.5*chi1(1)/avg.dTdz(iwind)^2;
                avg.kstart(iwind) = stats.k_start;
                avg.kstop(iwind)  = stats.k_stop;
                
                avg.n_iter(iwind) = length(chi1);
                
                if savespec==1
                    % 02/17/16 - AP - save spectra
                    
                    fspec=freq;
                    tpspec(iwind,:) = tp_power;                    
                    % observed wavenumber
                    ks(iwind,:) = k;                    
                    % observed spectra
                    kspec(iwind,:) = spec;                    
                    % best-fit theoreticl spectra
                    kkspec(iwind,:) = speck;                    
                    if ~isnan(kk)
                        % theoretical fit wavenumber
                        kks(iwind,:) = kk;
                    end
                end
                
            end % if T1Pvar>threshold
            
        end % windows
        
        
        if makeplots==1
            %~ plot summary figure
            ax = CTD_chipod_profile_summary(avg,cal,TP);
            axes(ax(1))
            title(Flist(icast).name,'interpreter','none')
            %~~~
        end
        
        if savespec==1
            %~
            avg.tpspec = tpspec;
            avg.kspec  = kspec;
            avg.kkspec = kkspec;
            avg.ks  = ks;
            avg.kks = kks;
            avg.fspec = fspec;
        end
        
        avg.lat = ctd.lat;
        avg.z   = sw_dpth(avg.P,avg.lat);
        avg.MakeInfo = ['Made ' datestr(now) ' w/ ComputeChi_Chameleon_tiwe.m'];
        
        % save results
        save( fullfile( datdirsave,[project_short '_' sprintf('%04d',icast) '_avg.mat']),'avg')
        
    catch
        disp(['error on profile ' num2str(icast)])
    end
    
end % icast

delete(hb)
%%