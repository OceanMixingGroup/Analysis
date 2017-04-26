function [chipod, cham] =Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compile data from binned chipod method and chameleon for specified
% profiles, averaged in bins of size dz. For eq08 or eq14.
%
% Similar to Get_binned_profiles, but returns matrix of profiles instead of
% one vector. (each profile is binned, but profiles are not averaged)
%
% - log10(chamleon epsilon) < -8.5 are discarded
%
%
% INPUT
% path_chipod_bin
% path_cham_avg
% dz
% Params
%
% OUTPUT
% - chipod : structure w/ binned profiles
% - cham   : ""
%
%-----------------
% 4/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Make empty arrays

empty_array = nan * ones(length([zmin:dz:zmax]),length(cnums_to_get) );

eps_chi  = empty_array ;
eps_cham = empty_array ;

chi_chi  = empty_array;
chi_cham = empty_array;

N2_chi  = empty_array;
N2_cham = empty_array;

Tz_chi  = empty_array;
Tz_cham = empty_array;

%P_cham_avg = empty_array;
%P_chi_avg  = empty_array;

hb = waitbar(0,['getting binned profiles for ' project_short])
for ic = 1:length(cnums_to_get)
    
    waitbar(ic/length(cnums_to_get),hb)
    
    clear cnum
    cnum = cnums_to_get(ic);
    
    clear avg ch chb
    
    try
        
        % regular chi-pod method on binned data
        clear avg
        if strcmp(project_short,'eq14')
            load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[upper(project_short) '_' sprintf('%04d',cnum) '_avg.mat']))
        else
            load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
        end
        chb = avg ; clear avg
        
        % chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        %% discard chameleon epsilons below noise floor
        
        clear ib
        ib = find( log10(avg.EPSILON)<-8.5 );
        avg.EPSILON(ib) = nan ;
        
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

chipod.P  = [zmin:dz:zmax]';
chipod.cnum = cnums_to_get ;
cham.P    = [zmin:dz:zmax]';
cham.cnum = cnums_to_get   ;

%%