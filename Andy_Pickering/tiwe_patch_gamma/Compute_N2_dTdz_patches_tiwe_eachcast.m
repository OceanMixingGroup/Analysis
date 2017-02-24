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

save_dir_patch='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/patches/'
addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

hb=waitbar(0)
warning off

for cnum=3818:4000
    waitbar(cnum/4000,hb)
    try
        
        % load patch data for this profile
        clear patch_data patches
        load(fullfile(save_dir_patch,['tiwe_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']))
        
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
        cham_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/cal';
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
        
        for ip=1:Npatches
            
            clear pstarts pstops
            clear iz t_ot s_ot p_ot ptmp_ot sgth_ot
            
            % get indices of data in patch
            iz=isin(p,[ patches.p1(ip) patches.p2(ip) ]);
            
            if length(iz)>10
                
                t_ot=t(iz);
                s_ot=s(iz);
                p_ot=p(iz);
                ptmp_ot=ptmp(iz);
                sgth_ot=sgth(iz);
                
                clear t_sort I
                [t_sort , I]=sort(t_ot,1,'descend');
                
                % sort potential temp
                clear t_pot t_pot_sort
                [ptmp_sort , Iptmp]=sort(ptmp_ot,1,'descend');
                
                clear DT dz dTdz
                dT=nanmax(ptmp_ot)-nanmin(ptmp_ot);
                dz=nanmax(p_ot)-nanmin(p_ot);
                dTdz=-dT/dz;
                
                % fit a line
                P=polyfit(p_ot,ptmp_sort,1);
                
                % save results
                patches.dtdz_range(ip)=dTdz;
                
                patches.dtdz_line(ip)=-P(1);
                
                %~~ 'bulk gradient' method from Smyth et al 2001
                % essentially = rms T (btw sorted/unsorted) /  thorpe scale ?
                %    t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
                t_rms= sqrt( nanmean(( ptmp_ot - ptmp_sort ).^2) );
                patches.dtdz_bulk(ip)=  t_rms / patches.Lt(ip) ;
                
                %~~ Now do similar for density / N^2
                
                [sgth_sort , I]=sort(sgth_ot,1,'ascend');
                
                % try the range/dz method
                drho=nanmax(sgth_ot)-nanmin(sgth_ot);
                dz=nanmax(p_ot)-nanmin(p_ot);
                n2_1=9.81/nanmean(sgth_ot)*drho/dz;
                
                patches.n2_range(ip)=n2_1;
                
                % fit a line to sgth
                clear P1
                P1=polyfit(p_ot,sgth_sort,1);
                
                % calculate N^2 from this fit
                clear drhodz n2_2 drho dz n2_3
                drhodz=-P1(1);
                patches.drhodz_line(ip)=drhodz;
                
                n2_2=-9.81/nanmean(sgth)*drhodz;
                patches.n2_line(ip)=n2_2;
                
                % Smyth eta al says associated N^2=g*bulk gradient (assuming
                % density controlled by temperature??)
                % I think missing divide by rho_0 also?
                %patches.n2_bulk(ip)= 9.81 / nanmean(sgth) * t_rms / patches.Lt(ip)    ;
                alpha = -0.2607; % from fit of sgth vs theta
                patches.n2_bulk(ip)= -9.81 / nanmean(sgth) * alpha * t_rms / patches.Lt(ip)    ;
                
                % try computing drho/dz using 'bulk' method (instead of
                % assuming rho only depends on T
                clear rho_rms
                rho_rms = sqrt( nanmean(( sgth_ot - sgth_sort ).^2) );
                patches.drhodz_bulk(ip) = - rho_rms / patches.Lt(ip) ;
                patches.n2_bulk_2(ip) = - 9.81 / nanmean(sgth) * patches.drhodz_bulk(ip);
                
                % compute N^2 w/ sw_bfrq
                clear n2
                n2=sw_bfrq(s_ot(I),t_ot(I),p_ot,0.5);
                patches.n4(ip)=nanmean(n2);
                
            end
            
        end % ip
        
        % save profile
       save(fullfile(save_dir_patch,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']), 'patches' ) 
    end % try
end % cnum

delete(hb)
warning on
%%