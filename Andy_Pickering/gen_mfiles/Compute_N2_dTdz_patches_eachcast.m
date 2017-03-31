function [] = Compute_N2_dTdz_patches_eachcast(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_patches_eachcast.m
%
%
% Compute N2 and dT/dz for overturns in chameleon patches using a few
% different methods.
%
% Uses patches found in Find_Patches_*project_name*_Raw.m
%
% OUTPUT:
% 'patches' structure
%
%
%
%--------------
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)];

% set paths
eval([project_name '_patches_paths'])

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

hb=waitbar(0,['Computing N2, dTdz for ' project_name ' patches'])
warning off
ic=0;
for cnum=cnums_to_do
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    try
        
        % load raw patch data for this profile
        clear patch_data patches
        if merge_patches==1
            load(fullfile(save_dir_patch,ot_dir,'raw_merge',[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']))
        else
            load(fullfile(save_dir_patch,ot_dir,'raw',[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']))
        end
        
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
        patches.n2_line=EmpVec;
        patches.n2_bulk=EmpVec;
        patches.n4=EmpVec;
        
        % density gradients
        patches.drhodz_bulk=EmpVec;
        patches.drhodz_line=EmpVec;
        
        % Different methods of computing dT/dz
        %patches.dtdz_range=EmpVec;
        patches.dtdz_line=EmpVec;
        patches.dtdz_bulk=EmpVec;
        
        patches.chi = EmpVec;
        patches.eps = EmpVec;
        
        patches.gam_line=EmpVec;
        patches.gam_bulk=EmpVec;
        patches.gam4 = EmpVec;
        
        % load raw chameleon cast
        clear cal cal2 head
        
        % Load the data for this cast
        eval(['cal = load_cal_' project_name '(cnum) ;'])
        
        % compute pot. temp, pot. density etc.
        clear s t p lat ptmp sgth
        s = cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
        %s = smooth(s,20);
        t = cal.T1 (1:end-1);
        p = cal.P  (1:end-1);
        ptmp=sw_ptmp(s,t,p,0);
        sgth=sw_pden(s,t,p,0);
        
        for ip=1:Npatches
            clear out
            out=compute_Tz_N2_for_patch(patches.p1(ip), patches.p2(ip) ,p...
                ,t ,s ,ptmp ,sgth , patches.Lt(ip) ) ;
            
            patches.dtdz_line(ip) = out.dtdz_line ;
            patches.dtdz_bulk(ip) = out.dtdz_bulk ;
            
            patches.n2_line(ip) = out.n2_line ;
            patches.n2_bulk(ip) = out.n2_bulk ;
            patches.n4(ip) = out.n4 ;

        end % ip
        
        % save profile
        if merge_patches==1
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
        else
            save(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' )
        end
        
    catch
        disp(['error on ' num2str(cnum)])
    end % try
end % cnum

delete(hb)
warning on
%%