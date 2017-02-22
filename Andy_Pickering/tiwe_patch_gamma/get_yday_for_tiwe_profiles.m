%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% get_yday_for_tiwe_profiles.m
%
% Get yearday for each tiwe cast, from the head.startime field
%
%------------
% 2/22/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

addpath /Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/mfiles
%addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
%addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
%addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

path_raw='/Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/cham/tw/';
%path_save='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/cal/';
%ChkMkDir(path_save)

Flist=dir( fullfile(path_raw, '*tw91*'))

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'tw91';
q.series={'fallspd','t1','t2','t','c','s','theta','sigma','epsilon1','epsilon2','chi'...
    'az2','ax_tilt','ay_tilt'};
warning off

cnum=[];
yday=[];

for cast=1:4000%1394%[7:3918]%[858:1219,2123:2590]%
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
            %cali_tw91;
            clear yd
            yd=str2num(head.starttime(end-5:end));
            
            cnum=[cnum(:) ; cast];
            yday=[yday(:) ; yd];
            
        end % try
    end % file exists     
end % cast

%%

save('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_cnum_yday.mat','cnum','yday')

%%

figure(1);clf
plot(cnum,yday)
xlabel('cnum')
ylabel('yday')

%%
