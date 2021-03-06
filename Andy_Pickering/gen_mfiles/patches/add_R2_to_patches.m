function []=add_R2_to_patches(project_name,patch_size_min,...
    usetemp,merge_patches,min_sep,cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Add R^2 for T-S fit in each patch.
%
%---------
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)];

% set paths
eval([project_name '_patches_paths'])

% load combined patch data for all profiles
eval(['patches = load_' project_name '_patches_comb(patch_size_min, usetemp, merge_patches, min_sep) ;'])
%%
addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

% add empty array for R^2, and N^2 for good patches
patches.R2 = nan*ones(size(patches.cnum));
patches.n2_line_fit = nan*ones(size(patches.cnum));
patches.gam_line_fit = nan*ones(size(patches.cnum));

hb=waitbar(0,['adding R^2 to patches for ' project_name])
ic=0;
for cnum=cnums_to_do 
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    try
        
        % find patches for this profile
        clear igc Npatches EmpVec
        
        igc = find(patches.cnum==cnum);
        
        Npatches = length(igc);
        
        % load raw chameleon cast
        clear cal cal2 head
        eval(['cal = load_cal_' project_name '(cnum);']);
        cnum_loaded = cnum;
        
        for ip=1:Npatches
            
            clear t sraw idz p P S
            idz  = isin(cal.P,[patches.p1(igc(ip)) patches.p2(igc(ip))]);
            t    = cal.T1(idz) ;
            sraw = cal.SAL(idz);
            p    = cal.P(idz) ;
            [P,S]= polyfit(t, sraw, 1);
            
            makeplots = 0 ;
            if makeplots==1
                figure(1);clf
                
                subplot(121)
                plot(t,sraw)
                axis ij
                grid on
                xlabel('T')
                xlabel('P')
                
                subplot(122)
                plot(t,sraw,'o')
                grid on
                hold on
                plot(t,polyval(P,t))
                xlabel('T')
                ylabel('S')
                
            end
            
            % compute R^2 for fit
            R2 = compute_R2(t,sraw,P) ;
            
            patches.R2(igc(ip))=R2;
            
            %use fit to compute N^2
            
            sfit = polyval(P,t);
            ptmp=sw_ptmp(sfit,t,p,0);
            sgth=sw_pden(sfit,t,p,0);
            
            % sorth pot. dens.
            clear sgth_sort I
            [sgth_sort , I]=sort(sgth,1,'ascend');
            
            %~~ compute drho/dz and N^2
            % fit a line to sgth
            clear P1
            P1 = polyfit(p,sgth_sort,1);
            
            % calculate N^2 from this fit
            clear drhodz
            drhodz = -P1(1);
            patches.n2_line_fit(igc(ip)) = -9.81/nanmean(sgth)*drhodz;
            
        end % each patch
    catch
        disp(['error on cnum ' num2str(cnum)])
    end % try
    
end % cnum
delete(hb)

%% compute gamma for these values also
%patches.gam_line_fit = nan*(ones(size(patches.gam_bin))) ;
patches.gam_line_fit = ComputeGamma(patches.n2_line_fit, patches.dtdz_line, patches.chi, patches.eps);


%% re-save combined structure

if merge_patches==1
    save( fullfile( analysis_dir, project_long, 'data',...
        [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100) '.mat']), 'patches' )
else
    save( fullfile( analysis_dir, project_long, 'data',...
        [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )
end

%%