function [chipod, cham] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Get data from all chipod & cham profiles, then average in 10m
% bins., instead of binning each profile and then averaging profiles.
%
% * chameleon epsilons below noise floor (log10<-8.5) are discarded
%
% * some screeining done to chiod epsilons
%
% * also attempts to do some screening of spikes in chipod data
%
% Compile data from binned chipod method and chameleon for specified
% profiles, averaged in bins of size dz. Then average profiles together.
% For eq08 and eq14.
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
% - Pmin            : All data where (P < Pmin) nan'd out
% - screen_chi      : Nan chipod (log) epsilons below -8.5 (noise floor of
% Chameleon)
% - screen_ml       : Nan out mixed layer depths that are convectively
% unstable (determined in Identify_mixedlayer_eq 08/14
%
% OUTPUT
%
% Returns single profiles of data that are averaged in depth bins across
% all the profiles.
%
%-----------------
% 4/12/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Make empty arrays

eps_chi  = [] ;
eps_cham = [] ;

chi_chi  = [];
chi_cham = [];

N2_chi  = [];
N2_cham = [];

Tz_chi  = [];
Tz_cham = [];

P_chi = [];
P_cham = [];

cnums = [];
%hb = waitbar(0,['averaging profiles'])
for cnum = cnums_to_get
    
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
        
        % load chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        % discard data in convectively unstable regions
        if screen_ml==1
            if strcmp(project_short,'eq14')
                avg = discard_convection_eq14_cham(avg,cnum);
            elseif strcmp(project_short,'eq08')
                avg = discard_convection_eq08_cham(avg,cnum);
            end
        end
        
        P_chi  = [ P_chi(:) ; chb.P(:) ] ;
        P_cham = [ P_cham(:) ; avg.P(:)] ;
        
        eps_cham = [eps_cham(:) ; avg.EPSILON(:) ];
        eps_chi  = [eps_chi(:)  ; chb.eps1(:) ];
        
        chi_cham = [chi_cham(:) ; avg.CHI(:) ];
        chi_chi  = [chi_chi(:)  ; chb.chi1(:) ];
        
        N2_cham = [N2_cham(:) ; avg.N2(:) ];
        N2_chi  = [N2_chi(:)  ; chb.N2(:) ];
        
        Tz_cham = [Tz_cham(:) ; avg.DTDZ(:) ];
        Tz_chi  = [Tz_chi(:)  ; chb.dTdz(:) ];
        
        cnums = [ cnums cnum];
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum
%delete(hb)

% Nan out values shallower than Pmin

clear ib
ib = find(P_cham<Pmin);
chi_cham(ib) = nan;
eps_cham(ib) = nan;

clear ib
ib = find(P_chi<Pmin);
chi_chi(ib) = nan;
eps_chi(ib) = nan;


% discard chameleon epsilons below noise floor
clear ib
ib = find(log10(eps_cham)<-8.5);
eps_cham(ib) = nan ;

ib = find(log10(eps_cham)>-5);
eps_cham(ib) = nan ;


% discard chipod epsilons below 8.5 (same as chanmeleon)
if screen_chi==1
    
    clear ib
    ib = find(log10(eps_chi)<-8.5);
    chi_chi(ib) = nan;
    eps_chi(ib) = nan;
    
    clear ib
    ib = find(log10(eps_chi)>-5);
    chi_chi(ib) = nan;
    eps_chi(ib) = nan;
    
end

%% now bin-average profiles together

[cham.eps z1 Nobs] = binprofile(eps_cham, P_cham, 0, dz, 200,1);
[cham.chi z1 Nobs] = binprofile(chi_cham, P_cham, 0, dz, 200,1);
[cham.N2  z1 Nobs] = binprofile(N2_cham,  P_cham, 0, dz, 200,1);
[cham.Tz  z1 Nobs] = binprofile(Tz_cham,  P_cham, 0, dz, 200,1);

cham.P = z1 ;
clear z1

[chipod.eps z1 Nobs] = binprofile(eps_chi, P_chi, 0, dz, 200,1);
[chipod.chi z1 Nobs] = binprofile(chi_chi, P_chi, 0, dz, 200,1);
[chipod.N2  z1 Nobs] = binprofile(N2_chi,  P_chi, 0, dz, 200,1);
[chipod.Tz  z1 Nobs] = binprofile(Tz_chi,  P_chi, 0, dz, 200,1);

chipod.P = z1;

cham.cnum = cnums ;
chipod.cnum = cnums;

%%