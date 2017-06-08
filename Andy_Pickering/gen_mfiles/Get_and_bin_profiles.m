function [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax,Pmin,screen_chi,screen_ml)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compile data from binned chipod method (applied to Chameleon thermistor
% profiles and chameleon for specified profiles, averaged in bins of size
% dz. For eq08 or eq14 data.
%
% Similar to Get_binned_profiles, but returns matrix of profiles instead of
% one vector. (each profile is binned, but profiles are not averaged)
%
% - log10(chamleon epsilon) < -8.5 are discarded
% 
%
%
% INPUT
% - path_chipod_bin : Set in eq08_patches_paths or eq14_patches_paths
% - path_cham_avg   : Set in eq08_patches_paths or eq14_patches_paths
% - dz              : bin size to average in
% - Params          : Params for chipod method 
% - cnums_to_get    : cast numbers to retrieve data for
% - project_short   : Project name (Set in eq08_patches_paths or
% eq14_patches_paths)
% - zmin            : Min P for binning profiles
% - zmax            : Max P for binning profiles
% - Pmin            : All data where (P < Pmin) nan'd out
% - screen_chi      : Nan chipod (log) epsilons below -8.5 (noise floor of
% Chameleon)
% - screen_ml       : Nan out mixed layer depths that are convectively
% unstable (determined in Identify_mixedlayer_eq 08/14
%
% OUTPUT
% - chipod : structure w/ binned profiles (matrix where columsn are
% profiles)
% - cham   : ""
%   - chipod,cham have fields: eps,chi,N2,dTdz
%
% DEPENDS ON
% - load_chipod_avg.m
% - discard_convection_eq14_cham.m
% - discard_convection_eq14_chi.m
% - binprofile.m
%
%-----------------
% 4/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Make empty arrays for results

empty_array = nan * ones(length([zmin:dz:zmax]),length(cnums_to_get) );

eps_chi  = empty_array ;
eps_cham = empty_array ;

chi_chi  = empty_array;
chi_cham = empty_array;

N2_chi  = empty_array;
N2_cham = empty_array;

Tz_chi  = empty_array;
Tz_cham = empty_array;

hb = waitbar(0,['getting binned profiles for ' project_short]);

for ic = 1:length(cnums_to_get)
    
    waitbar(ic/length(cnums_to_get),hb)
    
    clear cnum
    cnum = cnums_to_get(ic);
    
    clear avg ch chb
    
    try
               
        % load chipod-method profile
        chb = load_chipod_avg(path_chipod_bin,project_short,Params,cnum) ;
        
        % discard data in convectively unstable regions
        if screen_ml==1
            if strcmp(project_short,'eq14')
                chb = discard_convection_eq14_chi(chb,cnum);
            elseif strcmp(project_short,'eq08')
                chb = discard_convection_eq08_chi(chb,cnum);
            end
        end
        
        % NaN data shallower than Pmin
        izb = find(chb.P<Pmin);
        chb.eps1(izb) = nan;
        chb.chi1(izb) = nan;
        
        % load chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        % NaN data shallower than Pmin
        izb = find(avg.P<Pmin);
        avg.EPSILON(izb) = nan;
        avg.CHI(izb)     = nan;
        
        % discard data in convectively unstable regions
        if screen_ml==1            
            if strcmp(project_short,'eq14')
                avg = discard_convection_eq14_cham(avg,cnum);
            elseif strcmp(project_short,'eq08')
                avg = discard_convection_eq08_cham(avg,cnum);
            end
        end
        
        % discard chameleon epsilons below noise floor        
        clear ib
        ib = find( log10(avg.EPSILON)<-8.5 );
        avg.EPSILON(ib) = nan ;
        %avg.EPSILON(ib) = 1e-12 ;
        
        if screen_chi==1
            clear ib
            ib = find( log10(chb.eps1)<-8.5);
            chb.eps1(ib) = nan ;
            ib = find( log10(chb.eps1)>-5);
            chb.eps1(ib) = nan ;
            % chb.eps1(ib) = 1e-12 ;
        end
        
        [eps_cham(:,ic), ~ , ~] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,0);
        [eps_chi(:,ic) , ~ , ~] = binprofile(chb.eps1   , chb.P, 0, dz, 200,0);
        
        [chi_cham(:,ic), ~ , ~] = binprofile(avg.CHI    , avg.P, 0, dz, 200,0);
        [chi_chi(:,ic) , ~ , ~] = binprofile(chb.chi1   , chb.P, 0, dz, 200,0);
        
        [N2_cham(:,ic), ~ , ~] = binprofile(avg.N2   , avg.P, 0, dz, 200,0);
        [N2_chi(:,ic) , ~ , ~] = binprofile(chb.N2   , chb.P, 0, dz, 200,0);
        
        [Tz_cham(:,ic) ,~ , ~] = binprofile(avg.DTDZ , avg.P, 0, dz, 200,0);
        [Tz_chi(:,ic)  ,~ , ~] = binprofile(chb.dTdz   , chb.P, 0, dz, 200,0);
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum
delete(hb)

chipod = struct('eps',eps_chi, 'chi',chi_chi ,'N2',N2_chi ,'Tz',Tz_chi);
cham   = struct('eps',eps_cham,'chi',chi_cham,'N2',N2_cham,'Tz',Tz_cham);

chipod.P    = [zmin:dz:zmax]';
chipod.cnum = cnums_to_get ;
cham.P      = [zmin:dz:zmax]';
cham.cnum   = cnums_to_get   ;

%%