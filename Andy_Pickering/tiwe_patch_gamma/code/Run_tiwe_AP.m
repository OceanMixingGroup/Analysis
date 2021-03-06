%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Run_tiwe_AP.m
%
% Modified from run_tw91.m
%
%
% Dependencies:
%   - raw_load.m
%   - cali_tw91.m
%   - average_data_gen1.m
%
%--------------
% 2/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all ; clc

addpath /Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/mfiles
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate

path_raw  = '/Users/Andy/Dropbox/AP_Share_With_JN/date_from_jim/Tiwe91/cham/tw/' ;

tiwe_patches_paths
%%
ChkMkDir(path_cham_avg)

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'tw91';
q.series={'fallspd','t1','t2','t','c','s','theta','sigma','epsilon1','epsilon2','chi'...
    'az2','ax_tilt','ay_tilt'};
warning off

hb = waitbar(0,'processing tiwe chameleon files')

for cast = 1992:4000
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
            nfft=256;
            warning off
            if bad~=1
                avg=average_data_gen1(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
                
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
                
                % calc dynamic height
                head=calc_dynamic_z(avg,head);
                
                % add N2 and dT/dz to 'avg'
                avg.N2  = sw_bfrq(avg.S,avg.T,avg.P,0);
                avg.N2  = [avg.N2(:) ; nan]';
                avg.DTDZ= diffs(avg.T)./diffs(avg.P);
                avg.DTDZ= avg.DTDZ(:)';
                
                % create a seperate .mat data file containing 1m binned data and header
                temp=num2str(q.script.num+10000);
                fn=[q.script.prefix '_' temp(2:5) '_avg'];
                head.p_max=max(cal.P);
                %
                %eval(['save ' path_save fn ' avg head']);
                save(fullfile(path_cham_avg,fn),'avg','head')
            end % if not bad
        catch
            disp('error')
        end % try
    end % if file exists
end % cast

delete(hb)

%sum_tw91

%%
