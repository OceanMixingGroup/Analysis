%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakeCasts_eq08.m
%
% Does some of the pre-processing that was being repeated in
% ComputeChi_Chameleon_Eq08.m, and saves files that are now loaded by that
% script
%
% Formerly part of ComputeChi_Chameleon_Eq08.m
%
%----------------
% 5/15/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

makeplots = 0
savespec  = 1  % option to save wavenumber spectra

% Params for chipod calculations
Params.z_smooth = 1  ;  % distance to smooth N^2 and dTdz over

% specify method of computing N^2 and dT/dz
whN2dTdz  = 'regular'

% Add all the paths we need from mixing software
mixpath = '/Users/Andy/Cruises_Research/mixingsoftware/' ;
addpath(fullfile(mixpath,'seawater'))
addpath(fullfile(mixpath,'general'))
addpath(fullfile(mixpath,'marlcham'))
addpath(fullfile(mixpath,'CTD_Chipod','mfiles')) ;
addpath(fullfile(mixpath,'chipod','compute_chi')); % get_chipod_chi.md

eq08_patches_paths


%%
disp(['Data will be saved to ' path_chipod_bin])

% check if directory exists, make new if not
ChkMkDir(path_chipod_bin)
ChkMkDir(fullfile(path_chipod_bin,'cal',['zsmooth_' num2str(Params.z_smooth)]))
% initialize a waitbar
hb=waitbar(0);

%%
tstart=tic;
% loop through each cast

cnums_to_do = 193:2700 ;

warning off

for icast = 1:length(cnums_to_do)
    
    clear cal head ctd
    
    % update waitbar
    waitbar(icast/length(cnums_to_do),hb,['time elapsed=' num2str(round(toc(tstart)/60)) 'mins'])
    
    clear cnum
    cnum = cnums_to_do(icast)
    
    try
        
        close all
        clear cal head
        clear tpspec kspec kkspec fspec kks ks
        
        % Load the data for this cast
        load(fullfile(save_dir_cal,['eq08_' sprintf('%04d',cnum) '_avg.mat']))
        clear cal
        cal     = cal2;
        cal.T1  = cal.T;
        cal.SAL = cal.S ;
        cal.FALLSPD = cal.fspd;
        
        % Average temp and sal in 1m bins like we normally do for CTD data
        clear zmin dz zmax tbin zbin sbin
        zmin = nanmin(cal.P);
        dz   = 1;
        zmax = nanmax(cal.P);
        minobs = 2;
        [tbin zbin Nobs] = binprofile(cal.T1 ,cal.P, zmin, dz, zmax, minobs);
        [sbin zbin Nobs] = binprofile(cal.SAL,cal.P, zmin, dz, zmax, minobs);
        clear zmin dz zmax minobs
        
        % Put in 'ctd' structure
        clear ctd z_smooth
        ctd    = struct();
        ctd.t1 = tbin;
        ctd.s1 = sbin;
        ctd.p  = zbin;
        
        % add in lat and lon (from sum_eq14.m)
        clear idot lat1 lat2
        idot = strfind(head.lon.start,'.');
        lon1 = str2num(head.lon.start(1:idot-3));
        lon2 = str2num(head.lon.start(idot-2:end))/60;
        ctd.lon = nanmean([lon1 lon2]);
        
        clear idot lat1 lat2
        idot = strfind(head.lat.start,'.');
        lat1 = str2num(head.lat.start(1:idot-3));
        lat2 = str2num(head.lat.start(idot-2:end))/60;
        ctd.lat = nanmean([lat1 lat2]);
        %ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        
        ctd = Compute_N2_dTdz_forChi(ctd,Params.z_smooth);
        
        % save data here for loading during later chipod processing
        save( fullfile(path_chipod_bin,'cal',['zsmooth_' num2str(Params.z_smooth)],['eq08_' num2str(cnum) '_cal.mat']),'cal','ctd','head' )
        
    catch
        
    end % try
    
    
end % icast

warning on
%%