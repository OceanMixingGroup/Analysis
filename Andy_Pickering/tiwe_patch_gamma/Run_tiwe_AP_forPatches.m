%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Run_tiwe_AP_forPatches.m
%
% Run tiwe Chameleon processing for patches only. Patches are idenfitifed
% in FindPatches_tiwe_Raw.m, which must be run first.
% Compute_N2_dTdz_patches_tiwe.m must then be run also.
%
% Modified from Run_tiwe_AP.m
% Modified from run_tw91.m
%
% Dependencies:
%   - raw_load.m
%   - cali_tw91.m
%   - average_data_PATCH_AP.m
%
%--------------
% 2/16/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear all ; close all

patch_size_min = 0.15 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

addpath /Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/mfiles
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/ % average_data_PATCH_AP
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/mfiles/ % calc_chi_AP.m

path_raw='/Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/cham/tw/';

% set paths
tiwe_patches_paths

path_save = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
ChkMkDir(path_save)

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'tw91';
q.series={'fallspd','t1','t2','t','c','s','theta','sigma','epsilon1','epsilon2','chi'...
    'az2','ax_tilt','ay_tilt'};
warning off

hb=waitbar(0)

for cast=1:4000
    
    waitbar(cast/4000,hb)
    
    % bad files: 144
    %disp(cast);
    q.script.num=cast;
    temp1=q;
    clear global head data cal q
    global data head cal q
    q=temp1;
    raw_name=[q.script.pathname  q.script.prefix sprintf('%5.3f',q.script.num/1000)];
    if exist(raw_name,'file')==2
        try
            [data head]=raw_load(q);
            cali_tw91;
            
            % change S to SAL to conform w/ newer chameleon functions
            cal.SAL=cal.S;
            head.irep.SAL=head.irep.S;
            
            nfft=256;
            warning off
            if bad~=1
                
                % ~~ get patches for this profile
                clear patches
                % load the patches for this profile
                load(fullfile(save_dir_patch,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cast) '.mat']) )
                
                clear igp pstarts pstops
                pstarts = patches.p1 ;
                pstops  = patches.p2 ;
                
                % compute chi, eps etc. for patches
                %avg=average_data_gen_ct01a(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
                %avg=average_data_gen1(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
                avg=average_data_PATCH_AP(q.series,nfft,pstarts,pstops);
                
                % remove glitches
                % flag AZ vibrations
                %         idaz=find(avg.VARAZ>1.e-02);
                %         avg.EPSILON1(idaz)=NaN; avg.EPSILON2(idaz)=NaN; avg.EPSILON(idaz)=NaN;
                %         avg.WD(idaz)=NaN; avg.WD2(idaz)=NaN;
                tt=-5.3;
                temp2=find(log10(avg.AZ2)>tt);
                avg.EPSILON1(temp2)=NaN; avg.EPSILON2(temp2)=NaN; avg.EPSILON(temp2)=NaN;
                bad1=find(avg.EPSILON1>1e-4);
                avg.EPSILON1(bad1)=NaN;
                bad2=find(avg.EPSILON2>1e-4);
                avg.EPSILON2(bad2)=NaN;
                bad=find(avg.EPSILON>1e-4);
                avg.EPSILON(bad)=NaN;
                % flag surface
                idsur=find(avg.P<=5);
                avg.EPSILON1(idsur)=NaN; avg.EPSILON2(idsur)=NaN; avg.EPSILON(idsur)=NaN;
                %
                warning backtrace
                %
                % calc dynamic height
                %
                head=calc_dynamic_z(avg,head);
                
                
                % create a seperate .mat data file containing 1m binned data and header
                temp=num2str(q.script.num+10000);
                fn=[q.script.prefix temp(2:5) '_avg'];
                head.p_max=max(cal.P);
                %eval(['save ' path_save fn ' avg head']);
                save( fullfile(path_save,fn),'avg','head')
            end % if not bad
        end % try
    end % if file exists
end % cast
delete(hb)
%sum_tw91

%%
