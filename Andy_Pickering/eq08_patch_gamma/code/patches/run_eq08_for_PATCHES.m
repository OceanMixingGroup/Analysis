function [] = run_eq08_for_PATCHES(patch_size_min,usetemp,...
    merge_patches,min_sep,cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% run_eq08_for_PATCHES.m
%
% Script to run processing of EQ08 Chameleon data *for patches only*
%
% Modified from run_eq08.m, (by Sasha?).
%
% Dependencies:
% tag_file_eq08.m
% raw_load.m
% cali_eq08.m
%
%------------------
% 2/27/17 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%clear all

addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/general
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/
addpath /Volumes/'SP PHD U3'/NonBackup/EQ08/mfiles/

% folder with raw Chameleon data files
path_raw='/Volumes/SP PHD U3/NonBackup/EQ08/raw/'

eq08_patches_paths

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)]

% path to save processed files to
path_save = save_dir_avg_patch

% check if path_save exits, make new directory if not
ChkMkDir(path_save)

if merge_patches==1
    path_save = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_merged_minsep_' num2str(min_sep*100)])
else
    path_save = fullfile( save_dir_avg_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])
end
ChkMkDir(path_save)


% run script to tag bad casts and depth ranges before processing
tag_file_eq08

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'eq08';
series={'fallspd','t','tp','c','mht','s','theta','sigma','sigma_order','n2','epsilon1','epsilon2',...
    'chi','az2','dtdz','drhodz','varaz','varlt','varlt_theta','scat','ax_tilt','ay_tilt'};
warning off
ic=0;
hb=waitbar(0,'processing eq08 chameleon for patches')

for cast = cnums_to_do
    try
        %     [193:543,545:547,549,551:685,687:719,721:789,791,...
        %         793:798,800,802:988,990:1083,1085:1189,1191:1199,1201:1414,...
        %         1416:1420,1422:1500,1502,1504:1624,1631:1861,1863:1978,...
        %         1980:2124,2126:2140,2142:2298,2300:2668]
        %
        %
        
        %try
        %1524
        
        % bad, nonfixable files:
        % 79,792,799,801,1415,1421,1501,1503,1625,1862,2141,...
        %
        % up or surface profiles:
        % 544,548,550,686,720,989,1084,1190,1200,1979,2125,2299
        ic=ic+1;
        
        waitbar(ic/2624,hb)
        
        %disp(['working on cast ' num2str(cast)]);
        q.script.num=cast;
        q.series=series;
        temp1=q;
        clear global head data cal q
        global data head cal q
        q=temp1;
        
        % load raw Chameleon data
        [data head]=raw_load(q);
        
        % apply calibrations to data
        cali_eq08;
        warning off
        
        % * 3/15/16 - AP - skip this part, just want calibrated T' *
        %     average calibrated data into 1m bins
        %     nfft=128;
        nfft=256;
        
        % *** replacd with average_data_patches etc. ***
        % avg=average_data_gen1(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        % avg=average_data(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        
        %******
        % ~~ get patches for this profile
        clear patches
        % load the patches for this profile
        %    load(fullfile(save_dir_patch,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cast) '.mat']) )
        % load the patches for this profile
        clear patches
        if merge_patches==1
            load(fullfile(save_dir_patch, ot_dir, [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cast) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
        else
            load(fullfile(save_dir_patch,ot_dir, [project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cast) '.mat']) )
        end
        
        clear igp pstarts pstops
        pstarts = patches.p1 ;
        pstops  = patches.p2 ;
        
        % compute chi, eps etc. for patches
        avg=average_data_PATCH_AP(q.series,nfft,pstarts,pstops);
        
        %******
        
        
        %
        % remove glitches
        % flag AZ vibrations
        idaz=find(avg.VARAZ>1.e-02);
        avg.EPSILON1(idaz)=NaN; avg.EPSILON2(idaz)=NaN; avg.EPSILON(idaz)=NaN;
        %     avg.WD(idaz)=NaN; avg.WD2(idaz)=NaN;
        temp2=find(log10(avg.AZ2)>-4.5);
        avg.EPSILON1(temp2)=NaN; avg.EPSILON2(temp2)=NaN; avg.EPSILON(temp2)=NaN;
        bad1=find(avg.EPSILON1>1e-4);
        avg.EPSILON1(bad1)=NaN;
        bad2=find(avg.EPSILON2>1e-4);
        avg.EPSILON2(bad2)=NaN;
        bad=find(avg.EPSILON>1e-4);
        avg.EPSILON(bad)=NaN;
        %         avg.WD(temp2)=NaN; avg.WD2(temp2)=NaN;
        % flag surface
        idsur=find(avg.P<=5);
        avg.EPSILON1(idsur)=NaN; avg.EPSILON2(idsur)=NaN; avg.EPSILON(idsur)=NaN;
        if avg.EPSILON(end)>20*avg.EPSILON(end-1);
            avg.EPSILON(end)=NaN;avg.EPSILON(end)=NaN;avg.EPSILON(end)=NaN;
        end
        if q.script.num>=544 && q.script.num<=549
            avg.EPSILON=avg.EPSILON1;
        end
        %         avg.WD(idsur)=NaN; avg.WD2(idsur)=NaN;
        
        %             epsilon_glitch_factor=6;
        %             avg.EPSILON=(avg.EPSILON1+avg.EPSILON2)/2;
        %             % determine if EPSILON1>>>EPSILON2
        %             a=find(avg.EPSILON1>epsilon_glitch_factor*avg.EPSILON2 | isnan(avg.EPSILON1));
        %             avg.EPSILON(a)=avg.EPSILON2(a);
        %             % determine if EPSILON2>>>EPSILON1
        %             a=find(avg.EPSILON2>epsilon_glitch_factor*avg.EPSILON1 | isnan(avg.EPSILON2));
        %             avg.EPSILON(a)=avg.EPSILON1(a);
        
        
        warning backtrace
        
        % calc dynamic height
        head=calc_dynamic_z(avg,head);
        
        % save 'avg' data file containing 1m binned data and header
        temp=num2str(q.script.num+10000);
        fn=[q.script.prefix '_' temp(2:5) '.mat'];
        head.p_max=max(cal.P);
        %eval(['save ' path_save fn ' avg head']);
        %    save( fullfile(save_dir_avg_patch,fn),'avg','head')
        save( fullfile(path_save,fn),'avg','head')
        % end % try
        
    catch
        disp(['problem w/ cast ' num2str(cast)])
    end
end % cast #
delete(hb)

% sum_eq08_uppend
%sum_eq08

%%