function [chipod, cham] =Get_all_chipod_cham_data(path_chipod_bin,...
    path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml)
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
% - path_chipod_bin : Set in eq08_patches_paths or eq14_patches_paths
% - path_cham_avg   : Set in eq08_patches_paths or eq14_patches_paths
% - dz              : bin size to average in
% - Params          : Params for chipod method 
% - cnums_to_get    : cast numbers to retrieve data for
% - project_short   : Project name (Set in eq08_patches_paths or
% eq14_patches_paths)
% - Pmin            : All data where (P < Pmin) nan'd out
%
% OUTPUT
% - chipod
% - cham
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
        
        % load chipod-method profile
        chb = load_chipod_avg(path_chipod_bin,project_short,Params,cnum) ;
        
        % discard data in convectively unstable regions
        if screen_ml==1
            chb = discard_convection_eq14_chi(chb,cnum);
        end
        
        % Nan out values shallower than Pmin
        izb = find(chb.P<Pmin);
        chb.eps1(izb) = nan;
        chb.chi1(izb) = nan;
        
        
        if screen_chi==1
            clear ib
            ib = find( log10(chb.eps1)<-8.5 );
            chb.chi1(ib) = nan ;
            chb.eps1(ib) = nan ;
            
            clear ib
            ib = find( log10(chb.eps1)>-5 );
            chb.chi1(ib) = nan ;
            chb.eps1(ib) = nan ;
        end
                
        % load chamelon data (1m bins)
        load(fullfile( path_cham_avg, [project_short '_' sprintf('%04d',cnum) '_avg.mat']) )
        
        % discard data in convectively unstable regions
        if screen_ml==1
            avg = discard_convection_eq14_cham(avg,cnum);
        end
        
        
        %% discard chameleon epsilons below noise floor
        
        clear ib
        ib = find( log10(avg.EPSILON)<-8.5 );
        avg.CHI(ib) = nan ;
        avg.EPSILON(ib) = nan ;
        
        clear ib
        ib = find( log10(avg.EPSILON)>-5 );
        avg.CHI(ib) = nan ;
        avg.EPSILON(ib) = nan ;
        
        % Nan out values shallower than Pmin
        izb = find(avg.P<Pmin);
        avg.EPSILON(izb) = nan;
        avg.CHI(izb)     = nan;
                
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