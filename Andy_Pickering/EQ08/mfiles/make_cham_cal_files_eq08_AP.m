%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% make_cham_cal_files_eq08_AP.m
%
% Used to be part of run_eq08_AP
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

addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/general
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths

if use_15hz==1
    save_dir_cal = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc_AP_15hz/cal/';
elseif use_10hz==1
    save_dir_cal = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/EQ08/Data/cham_proc_AP_10hz/cal/';
end

ChkMkDir(save_dir_cal)

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
        
        %~~ check if TP,p,c,t sampled at same rate; if not, interpolate to same
        % length; I want all to be same length as TP
        %%
        irP   = head.irep.P ;
        irTP  = head.irep.TP;
        irT   = head.irep.T ;
        irS   = head.irep.S ;
        irFspd= head.irep.FALLSPD;
        %%
        if head.irep.P~=irTP
            clear P2
            % make new P vector length of TP (want to preserve TP higher
            % sampling rate)
            P2=nan*ones(length(cal.TP),1);
            P2( 1 : irTP/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            cal2.P=P2;
            cal2.TP=cal.TP;
        end
        
        if irT~=irTP
            clear P2
            P2=nan*ones(length(cal.T),1);
            P2( 1 : irT/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.T=interp1(P2(ig),cal.T(ig),cal2.P);
        end
        
        if head.irep.S~=irTP
            clear P2
            P2=nan*ones(length(cal.S),1);
            P2( 1 : irS/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.S=interp1(P2(ig),cal.S(ig),cal2.P);
        end
        
        if head.irep.FALLSPD~=irTP
            clear P2
            P2=nan*ones(length(cal.FALLSPD),1);
            P2( 1 : irFspd/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.fspd=interp1(P2(ig),cal.FALLSPD(ig),cal2.P);
            
        end
        %~~
        
        cal.MakeInfo  = ['Made by AP ' datestr(now) ' w/ make_cham_cal_files_eq08_AP.m'] ;
        cal2.MakeInfo = ['Made by AP ' datestr(now) ' w/ make_cham_cal_files_eq08_AP.m'] ;
        
        % save calibrated cast files here
        clear temp fn
        temp = num2str(q.script.num+10000);
        fn   = [q.script.prefix '_' temp(2:5) '_cal.mat'];
        save(fullfile(save_dir_cal,fn),'cal2','head')
        
    catch
        disp(['error on cast ' num2str(cast)])
    end % try
    
end % cast #
delete(hb)

%%