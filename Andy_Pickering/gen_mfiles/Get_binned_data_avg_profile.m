function [eps_cham_avg, chi_cham_avg, N2_cham_avg, Tz_cham_avg, eps_chi_avg, chi_chi_avg, N2_chi_avg, Tz_chi_avg] =...
    Get_binned_data_avg_profile(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
%
%
%-----------------
% 4/12/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%


% Make empty arrays

eps_chi  = [] ;
eps_cham = [] ;

chi_chi = [];
chi_cham = [];

N2_chi = [];
N2_cham = [];

Tz_chi = [];
Tz_cham = [];

%hb = waitbar(0,['averaging profiles'])
for cnum = cnums_to_get 
    
    clear avg ch chb
    
    try
        
        % regular chi-pod method on binned data
        clear avg
        if project_short=='eq14'
        load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[upper(project_short) '_' sprintf('%04d',cnum) '_avg.mat']))            
        else
        load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
        end
        chb = avg;clear avg
        
        % chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        clear bin1 bin2
        [bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
        
        eps_cham = [eps_cham    bin1(:) ];
        eps_chi  = [eps_chi     bin2(:) ];
        
        clear bin1 bin2
        [bin1 z1 Nobs] = binprofile(avg.CHI, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.chi1   , chb.P, 0, dz, 200,1);
        
        chi_cham = [chi_cham    bin1(:) ];
        chi_chi  = [chi_chi     bin2(:) ];
        
        clear bin1 bin2
        [bin1 z1 Nobs] = binprofile(avg.N2, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.N2   , chb.P, 0, dz, 200,1);
        
        N2_cham = [N2_cham    bin1(:) ];
        N2_chi  = [N2_chi     bin2(:) ];
        
        clear bin1 bin2
        [bin1 z1 Nobs] = binprofile(avg.DTDZ, avg.P, 0, dz, 200,1);
        [bin2 z2 Nobs] = binprofile(chb.dTdz   , chb.P, 0, dz, 200,1);
        
        Tz_cham = [Tz_cham    bin1(:) ];
        Tz_chi  = [Tz_chi     bin2(:) ];
        
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum
%delete(hb)
%% now average profiles together

eps_cham_avg = nanmean(eps_cham,2);
chi_cham_avg = nanmean(chi_cham,2);
N2_cham_avg  = nanmean(N2_cham,2) ;
Tz_cham_avg  = nanmean(Tz_cham,2) ;

eps_chi_avg = nanmean(eps_chi,2);
chi_chi_avg = nanmean(chi_chi,2);
N2_chi_avg  = nanmean(N2_chi,2) ;
Tz_chi_avg  = nanmean(Tz_chi,2) ;


%%