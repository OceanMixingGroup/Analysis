%~~~~~~~~~~~~~~~~~~~~~~~
%
% Process_tiwe_rawprofiles_AP.m
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

path_raw='/Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/cham/tw/';
path_save='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/cal/';
ChkMkDir(path_save)

Flist=dir( fullfile(path_raw, '*tw91*'))
%

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'tw91';
q.series={'fallspd','t1','t2','t','c','s','theta','sigma','epsilon1','epsilon2','chi'...
    'az2','ax_tilt','ay_tilt'};
warning off
for cast=2836:3711%1:4000%1394%[7:3918]%[858:1219,2123:2590]%
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
                    %cal.P=interp1(1:length(cal.P),cal.P,1:length(cal.TP));cal.P=cal.P(:);
                    P2=nan*ones(length(cal.T1P),1);
                    P2(1:head.irep.T1P:end)=cal.P;
                    P2=NANinterp(P2);
                    cal2.P=P2;
                    cal2.TP=cal.T1P;
                end
                
                if head.irep.T1~=irTP
                    %cal.T1=interp1(1:length(cal.T1),cal.T1,1:length(cal.TP));cal.T1=cal.T1(:);
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
                cal2.SAL_sm=smooth(cal2.SAL,20);
                
                %~~ save data
                temp=num2str(q.script.num+10000);
                fn=[q.script.prefix temp(2:5) '_raw'];
                head.p_max=max(cal.P);
                %disp('test')
                eval(['save ' path_save fn ' cal2 head']);
            end % if not bad
        end % try
    end % if file exists
end % cast
%sum_tw91

%%


%%