%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakeCasts_eq14.m
%
% Move repeated parts from ComputeChi_Chameleon_eq14.m here and save
% results for loading by that script.
%
%--------------
% 5/16/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Params for chipod calculations
Params.z_smooth = 10;    % distance to smooth N^2 and dTdz over
% Params.fmax     = 5;    % freq. to integrate dT/dz spectrum up to (ie noise level)
% Params.nfft     = 128;   % # points to use in computing spectra
% Params.TPthresh = 1e-6;  % minimum dT/dz variance
% Params.resp_corr= 0;     % correct TP spectra for freq response of thermistor?
% Params.fc       = 99 ;    % cutoff frequency for response correction
% Params.gamma    = 0.2 ;   % mixing efficiency

makeplots = 0
savespec  = 0  % option to save wavenumber spectra

% specify method of computing N^2 and dT/dz
whN2dTdz = 'regular'
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

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/code
eq14_patches_paths

% if Params.resp_corr==0
%     Params.fc=99;
% end

% Make directory to save processed casts in (name based on Params)

% if strcmp(whN2dTdz, 'regular')
%     datdirsave=fullfile(path_chipod_bin,...
%         ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft)]);
% else
%     datdirsave=fullfile(path_chipod_bin,...
%         ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft) '_' whN2dTdz]);
% end

datdirsave=fullfile(path_chipod_bin);

disp(['Data will be saved to ' datdirsave])

% check if directory exists, make new if not
ChkMkDir(datdirsave)

ChkMkDir(fullfile(datdirsave,'cal',['zsmooth_' num2str(Params.z_smooth)]))


% initialize a waitbar
hb=waitbar(0);

tstart=tic;
% loop through each cast

cnums_to_do = get_cham_cnums_eq14 ;
%cnums_to_do = 1000:2000 ;
%id = find(cnums_to_do>1650);
%cnums_to_do = cnums_to_do(id);

%%

for icast=1:length(cnums_to_do)
    
    cnum = cnums_to_do(icast);
    
    % update waitbar
    waitbar(icast/length(cnums_to_do),hb,['time elapsed=' num2str(round(toc(tstart)/60)) 'mins'])
    
    close all
    clear cal head
    clear tpspec kspec kkspec fspec kks ks
    
    try
        
        % Load the data for this cast
        load(fullfile(save_dir_cal,['EQ14_' sprintf('%04d',cnum) '.mat']))
        clear cal
        cal=cal2;
        
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
        clear idot lat1 lat2
        idot=strfind(head.lon.start,'.');
        lon1=str2num(head.lon.start(1:idot-3));
        lon2=str2num(head.lon.start(idot-2:end))/60;
        ctd.lon=nanmean([lon1 lon2]);
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        ctd.lat=nanmean([lat1 lat2]);
        % ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        
        % compute N^2 and dT/dz
        %        if strcmp(whN2dTdz, 'regular')
        ctd = Compute_N2_dTdz_forChi(ctd,Params.z_smooth);
        %         elseif strcmp(whN2dTdz, 'regular2')
        %             ctd = Compute_N2_dTdz_forChi_2(ctd);
        %         elseif strcmp(whN2dTdz, 'line')
        %             ctd = Compute_N2_dTdz_forChi_line(ctd,1);
        %         elseif strcmp(whN2dTdz, 'raw_line')
        %             % use raw (24hz) CTD data to comptue N^2 and Tz
        %             ctd = Compute_N2_dTdz_forChi_raw_line(ctd,cal,0.5);
        %         end
        %
        
        % save data here for loading during later chipod processing
        save( fullfile(datdirsave,'cal',['zsmooth_' num2str(Params.z_smooth)],['eq14_' num2str(cnum) '_cal.mat']),'cal','ctd','head' )
        
        
    catch
    end % try
    
end % icast



%%