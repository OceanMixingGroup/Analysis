%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% run_eq08_avg_AP.m
%
% Script to run processing of EQ08 Chameleon data.
%
% Modified from run_eq08.m, (by Sasha?).
%
% Dependencies:
% tag_file_eq08.m
% raw_load.m
% cali_eq08.m
% average_data_gen1.m
% average_data_gen1_AP_10hz.m
%
%------------------
% 03/15/16 - A.Pickering - apickering@coas.oregonstate.edu
% 03/18/16 - AP - Modifying to save N2 and dTdz also for use in chi calcs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

use_15hz = 0 % use fmax=15 instead
use_10hz = 1

% folder with raw Chameleon data files
path_cham_raw='/Volumes/SP PHD U3/NonBackup/EQ08/raw/'

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/general
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

eq08_patches_paths


if use_15hz==1
    path_cham_avg ='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc_AP_15hz/avg/';
elseif use_10hz==1
    path_cham_avg ='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc_AP_10hz/avg/';
end

ChkMkDir(path_cham_avg)

% run script to tag bad casts and depth ranges before processing
tag_file_eq08

global data head cal q
q.script.pathname =  path_cham_raw;
q.script.prefix = 'eq08';
series={'fallspd','t','tp','c','mht','s','theta','sigma','sigma_order','n2','epsilon1','epsilon2',...
    'chi','az2','dtdz','drhodz','varaz','varlt','varlt_theta','scat','ax_tilt','ay_tilt'};
warning off
ic=0;
hb=waitbar(0);

for cast=[193:543,545:547,549,551:685,687:719,721:789,791,...
        793:798,800,802:988,990:1083,1085:1189,1191:1199,1201:1414,...
        1416:1420,1422:1500,1502,1504:1624,1631:1861,1863:1978,...
        1980:2124,2126:2140,2142:2298,2300:2355,2357:2668]
    
    % bad, nonfixable files:
    % 79,792,799,801,1415,1421,1501,1503,1625,1862,2141,...
    %
    % up or surface profiles:
    % 544,548,550,686,720,989,1084,1190,1200,1979,2125,2299
    ic=ic+1;
    
    waitbar(ic/2624,hb)
    
    disp(['working on cast ' num2str(cast)]);
    
    try
        
        q.script.num = cast ;
        q.series = series ;
        temp1 = q ;
        clear global head data cal q
        global data head cal q
        q=temp1;
        
        % load raw Chameleon data
        [data head]=raw_load(q);
        
        % apply calibrations to data
        cali_eq08;
        warning off
        
        
        %     nfft=128;
        nfft=128;
        if use_15hz==1
            avg=average_data_gen1_AP_15hz(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        elseif use_10hz==1
            avg=average_data_gen1_AP_10hz(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        else
            avg=average_data_gen1(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        end
        %avg=average_data(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        %
        % remove glitches
        % flag AZ vibrations
        idaz = find(avg.VARAZ>1.e-02);
        avg.EPSILON1(idaz)=NaN; avg.EPSILON2(idaz)=NaN; avg.EPSILON(idaz)=NaN;
        %     avg.WD(idaz)=NaN; avg.WD2(idaz)=NaN;
        
        temp2 = find(log10(avg.AZ2)>-4.5);
        avg.EPSILON1(temp2) = NaN; avg.EPSILON2(temp2)=NaN; avg.EPSILON(temp2)=NaN;
        
        bad1 = find(avg.EPSILON1>1e-4);
        avg.EPSILON1(bad1)=NaN;
        
        bad2 = find(avg.EPSILON2>1e-4);
        avg.EPSILON2(bad2)=NaN;
        
        bad = find(avg.EPSILON>1e-4);
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
        
        %create a seperate .mat data file containing 1m binned data and header
        
        temp = num2str(q.script.num+10000);
        fn = [q.script.prefix '_' temp(2:5) '_avg.mat'];
        %fn=[q.script.prefix '_' temp(2:5) '_avg.mat'];
        head.p_max=max(cal.P);
        %    eval(['save ' path_cham_avg fn ' avg head']);
        save(fullfile(path_cham_avg,fn),'avg','head')
        
    catch
        disp(['error on cast ' num2str(cast)])
    end % try
    
end % cast #
delete(hb)

% sum_eq08_uppend
%sum_eq08

%%