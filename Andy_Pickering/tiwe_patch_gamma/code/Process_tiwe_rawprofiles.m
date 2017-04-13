%~~~~~~~~~~~~~~~~~~~~~~~
%
% Process_tiwe_rawprofiles_AP.m
%
% Process tiwe chameleon raw profiles into mat files ('cal') with
% calibrated t,s,TP etc.. These are *NOT* 1m averaged yet (I need these to
% apply the chi-pod method to the raw TP profiles.
%
% Processed files are saved in /save_dir_cal/, which is specified in
% tiwe_patches_paths.m
%
% Modified from Run_tiwe_AP.m
%
%------------------
% 2/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/mfiles
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

path_raw  = '/Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/cham/tw/';

tiwe_patches_paths
%%
ChkMkDir(save_dir_cal)
%%

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'tw91';
q.series={'fallspd','t1','t2','t','c','s','theta','sigma','epsilon1','epsilon2','chi'...
    'az2','ax_tilt','ay_tilt'};
warning off

hb = waitbar(0,'processing raw tiwe files')

for cast=3883:4000%
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
            
            if bad~=1
                
                % make a new modified structure that I will use for chipod method
                cal2=struct();
                
                %~~ check if TP,p,c,t sampled at same rate; if not, interpolate to same
                % length; I want all to be same length as TP
                irTP=head.irep.T1P;
                if head.irep.P~=irTP
                    clear P2
                    P2=nan*ones(length(cal.T1P),1);
                    P2(1:head.irep.T1P:end)=cal.P;
                    P2=NANinterp(P2);
                    cal2.P=P2;
                    cal2.TP=cal.T1P;
                end
                
                if head.irep.T1~=irTP                    
                    cal2.T1=interp1(cal.P,cal.T1,cal2.P);
                else
                    cal2.T1=cal.T1;
                end
                
                if head.irep.S~=irTP
                    %cal.SAL=interp1(1:length(cal.SAL),cal.SAL,1:length(cal.TP));cal.SAL=cal.SAL(:);
                    cal2.SAL=interp1(cal.P,cal.S,cal2.P);
                end
                
                if head.irep.FALLSPD~=irTP
                    % cal.FALLSPD=interp1(1:length(cal.FALLSPD),cal.FALLSPD,1:length(cal.TP));cal.FALLSPD=cal.FALLSPD(:);
                    cal2.FALLSPD=interp1(cal.P,cal.FALLSPD,cal2.P);
                end
                
                % salinity is spiky, so also save a smoothed profile for
                % future use (ie in identifying patches)
                %cal2.SAL_sm=smooth(cal2.SAL,20);
                
                %~~ save data
                temp=num2str(q.script.num+10000);
                fn=[q.script.prefix '_' temp(2:5) '_cal.mat'];
                head.p_max=max(cal.P);
                save(fullfile(save_dir_cal,fn),'cal2','head')
            end % if not bad
        catch
        end % try
    end % if file exists
end % cast
delete(hb)
%sum_tw91

%%


%%