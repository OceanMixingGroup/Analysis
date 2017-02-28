%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_patches_tiwe_eachcast.m
%
% Modified from Compute_N2_dTdz_patches_tiwe.m
% Does/saves each cast separately instead of all casts/patches. Switching
% to this because it kept crashing and I had to start over...
%
% Compute N2 and dT/dz for overturns in tiwe chameleon profiles using a few
% different methods.
%
% Uses patches found in Find_Patches_tiwe_Raw.m
%
% OUTPUT:
% 'patches' structure
%
% See also LookAtProfilesOT.m
%
% - Run FindPatches_tiwe_Raw.m 1st to identify patches
%
%--------------
% 2/22/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1

% set paths
tiwe_patches_paths

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

alpha = -0.2607; % from fit of sgth vs theta

hb=waitbar(0)
warning off

for cnum=1:4000
    waitbar(cnum/4000,hb)
    try
        
        % load patch data for this profile
        clear patch_data patches
        load(fullfile(save_dir_patch,[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']))
        
        patches = struct() ;
        patches.cnum = patch_data(:,1) ;
        patches.p1   = patch_data(:,2) ;
        patches.p2   = patch_data(:,3) ;
        patches.n2_ot= patch_data(:,5) ;
        patches.Lt   = patch_data(:,6) ;
        patches.yday = patch_data(:,7) ;
        
        % Make empty arrays for results
        Npatches=length(patches.p1);
        EmpVec=nan*ones(Npatches,1);
        
        % Different methods of computing N^2
        patches.n2_range=EmpVec;
        patches.n2_line=EmpVec;
        patches.n2_bulk=EmpVec;
        patches.n4=EmpVec;
        patches.n2_bulk_2=EmpVec;
        
        % density gradients
        patches.drhodz_bulk=EmpVec;
        patches.drhodz_line=EmpVec;
        
        % Different methods of computing dT/dz
        patches.dtdz_range=EmpVec;
        patches.dtdz_line=EmpVec;
        patches.dtdz_bulk=EmpVec;
        
        % load raw chameleon cast
        clear cal cal2 head
        %cham_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/cal';
        cham_dir = save_dir_cal
        load( fullfile( cham_dir, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )
        cal=cal2 ; clear cal2
        cnum_loaded = cnum;
        
        % compute pot. temp, pot. density etc.
        clear s t p lat ptmp sgth
        s = cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
        s = smooth(s,20);
        t = cal.T1 (1:end-1);
        p = cal.P  (1:end-1);
        ptmp=sw_ptmp(s,t,p,0);
        sgth=sw_pden(s,t,p,0);
        
        % get latitude for profile
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        lat=nanmean([lat1 lat2]);
        %
        for ip=1:Npatches
        
            clear out
            out=compute_Tz_N2_for_patch(patches.p1(ip), patches.p2(ip) ,p...
                ,t ,s ,ptmp ,sgth ,alpha , patches.Lt(ip) ) ;c
            
            patches.dtdz_range(ip) = out.dtdz_range ;
            patches.dtdz_line(ip) = out.dtdz_line ;
            patches.dtdz_bulk(ip) = out.dtdz_bulk ;
            
            patches.n2_range(ip) = out.n2_range ;
            patches.n2_line(ip) = out.n2_line ;
            patches.n2_bulk(ip) = out.n2_bulk ;
            patches.n4(ip) = out.n4 ;
            
        end % ip
        
        % save profile
        save(fullfile(save_dir_patch,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
    end % try
end % cnum

delete(hb)
warning on
%%