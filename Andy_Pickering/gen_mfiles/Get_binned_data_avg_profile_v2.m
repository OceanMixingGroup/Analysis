function [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg, P_chi, P_cham] =...
    Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin)
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
% path_chipod_bin
% path_cham_avg
% dz
% Params
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

%hb = waitbar(0,['averaging profiles'])
for cnum = cnums_to_get 
    
    clear avg ch chb
    
    try
        
        % regular chi-pod method on binned data
        clear avg
        if strcmp(project_short,'eq14')
        load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[upper(project_short) '_' sprintf('%04d',cnum) '_avg.mat']))            
        else
        load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
        end
        chb = avg;clear avg
        
        % chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
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
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum
%delete(hb)

%% Nan out values shallower than Pmin

clear ib
ib = find(P_cham<Pmin);
eps_cham(ib) = nan;

clear ib
ib = find(P_chi<Pmin);
eps_chi(ib) = nan;


%% discard chameleon epsilons below noise floor

clear ib
ib = find(log10(eps_cham)<-8.5);
eps_cham(ib) = nan ;

%% Nan out spikes in chi-pod Tz,epsilon
% 
% clear ib
% ib = find( medfilt1(Tz_chi,5) ./ Tz_chi  >2 ) ;
% eps_chi(ib)=nan;
% chi_chi(ib)=nan;
% Tz_chi(ib) =nan;
% N2_chi(ib) =nan;
% 
% clear ib
% ib = find(log10(eps_chi)>-4);
% eps_chi(ib) = nan;

%% now bin-average profiles together

[eps_cham_avg z1 Nobs] = binprofile(eps_cham, P_cham, 0, dz, 200,1);
[chi_cham_avg z1 Nobs] = binprofile(chi_cham, P_cham, 0, dz, 200,1);
[N2_cham_avg  z1 Nobs] = binprofile(N2_cham,  P_cham, 0, dz, 200,1);
[Tz_cham_avg  z1 Nobs] = binprofile(Tz_cham,  P_cham, 0, dz, 200,1);

P_cham = z1 ;
clear z1

[eps_chi_avg z1 Nobs] = binprofile(eps_chi, P_chi, 0, dz, 200,1);
[chi_chi_avg z1 Nobs] = binprofile(chi_chi, P_chi, 0, dz, 200,1);
[N2_chi_avg  z1 Nobs] = binprofile(N2_chi,  P_chi, 0, dz, 200,1);
[Tz_chi_avg  z1 Nobs] = binprofile(Tz_chi,  P_chi, 0, dz, 200,1);

P_chi = z1;
%%