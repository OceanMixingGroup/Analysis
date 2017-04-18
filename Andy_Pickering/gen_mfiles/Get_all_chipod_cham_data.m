function [chipod, cham] =Get_all_chipod_cham_data(path_chipod_bin,...
    path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Function to get all data from chipod and chameleon profiles (no binning)
% i.e. concats all profiles into one list for each variable
% 
% Works for eq08 or eq14.
%
% - NOTE log10(chamleon epsilon) < -8.5 are discarded
%
% INPUT
% -path_chipod_bin
% -path_cham_avg
% -dz
% -Params
%
% OUTPUT
%
%
%
%-----------------
% 4/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

P_cham = [];
P_chi  = [];

cnum_chi  = [] ;
cnum_cham = [];

for ic = 1:length(cnums_to_get)
    
    clear cnum
    cnum = cnums_to_get(ic);
    
    clear avg ch chb
    
    try
        
        % regular chi-pod method on binned data
        clear avg
        if project_short=='eq14'
            load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[upper(project_short) '_' sprintf('%04d',cnum) '_avg.mat']))
        else
            load( fullfile( path_chipod_bin, ['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
        end
        chb = avg ; clear avg
        
        % chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        %% discard chameleon epsilons below noise floor
        
        clear ib
        ib = find( log10(avg.EPSILON)<-8.5 );
        avg.EPSILON(ib) = nan ;
        
        P_chi  = [ P_chi(:)  ; chb.P(:) ] ;
        P_cham = [ P_cham(:) ; avg.P(:) ] ;
                
        eps_cham = [eps_cham(:) ; avg.EPSILON(:) ];
        eps_chi  = [eps_chi(:)  ; chb.eps1(:) ];
                 
        chi_cham = [chi_cham(:) ; avg.CHI(:) ];
        chi_chi  = [chi_chi(:)  ; chb.chi1(:) ];
                
        N2_cham = [N2_cham(:) ; avg.N2(:) ];
        N2_chi  = [N2_chi(:)  ; chb.N2(:) ];
        
        Tz_cham = [Tz_cham(:) ; avg.DTDZ(:) ];
        Tz_chi  = [Tz_chi(:)  ; chb.dTdz(:) ];

        cnum_cham = [ cnum_cham(:) ; cnum*ones(length(avg.P),1) ];
        cnum_chi  = [ cnum_chi(:)  ; cnum*ones(length(chb.P),1) ];
                                
    catch
        disp(['error on profile ' num2str(cnum) ])
    end % try
    
end % cnum

chipod = struct('eps',eps_chi ,'chi',chi_chi ,'N2',N2_chi ,'Tz',Tz_chi , 'P',P_chi, 'cnum', cnum_chi );
cham   = struct('eps',eps_cham,'chi',chi_cham,'N2',N2_cham,'Tz',Tz_cham, 'P', P_cham, 'cnum', cnum_cham );


%%