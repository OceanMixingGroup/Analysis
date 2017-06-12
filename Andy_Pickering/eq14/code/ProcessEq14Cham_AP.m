%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ProcessEq14Cham_AP.m
%
% Process raw Chameleon data files from Eq14 cruise. 
%
% My goal is to test the
% chipod methods by applying them to TP from Chameleon w/o using shear, and
% comparing to the 'True' chi/epsilon computed using the shear.
%
% This file was modified from Sally Warner's run_eq14.m (thanks Sally for
% commenting your code well!).
%
% OUTPUT
% 'cal','cal2' structures containing calibrated temp,sal profiles at high
% resolution. This script does NOT compute chi/epsilon or do any binning.
%
% DEPENDENCIES
% - raw_load_cham2
% - cali_eq14
%
%---------------------------
% 01/19/16 - A.Pickering - apickering@coas.oregonstate.edu
% 01/25/16 - AP - Fixed bug with interpolating P,T1,SAL to same length as
% TP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% add folder containing Chameleon processing routines
addpath /Users/Andy/Cruises_Research/mixingsoftware/Chameleon2/Version2015/

% this folder has calibrate functions
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/

% this folder has mfiles that Sally sent me with data, including script to
% do calibrations etc
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/mfiles/

% define important data and paths
cruise_id = 'EQ14';
depth_max = 250;
path_cham = '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/'
path_raw  = fullfile(path_cham,'raw')

% make directory for processed files
path_save='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/'
ChkMkDir(path_save)

%% define imporant variables
% define global variables which are called by many of the scripts
global data head cal q

% make structure q.script which will be used eventually to call the correct
% file by raw_load
q.script.pathname =  path_raw;
q.script.prefix = cruise_id;

% define the variables to be processed
% (note Pavan changed a number of the variable names, so they are different
% from the dynamo versions of the code)
% (note: neet to have T1, T2, SAL before epsilon and chi)
% (note: no need to list dT/dz or N2 here. They are calculated from 1m bin
% averaged data and therefore are not needed in this list)
series = {'fallspd','t1','t2','tp','cond','sal','theta','sigma',...
    'epsilon1','epsilon2',...
    'chi','az2','scat','ax_tilt','ay_tilt',...
    'varaz','varlt'};

warning off


%% PROCESS EACH CAST

% manually define the cast numbers to be processed
% (note: in the realtime code, the "bad" flag is used for casts that are no
% good for whatever reason. Here, you're just supposed to have a matrix so
% the bad files are skipped over. In other words, "bad" isn't implemented.)
for cast = 4%[4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
    %1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];
    
    try
        
        disp(['processing cast number: ' num2str(cast)]);
        q.script.num=cast;
        q.series=series;
        temp1=q;
        
        % weird cludgy fix for redefining q for each cast (need to delete the
        % other variables that are saved in q for previous casts)
        clear global head data cal q
        global data head cal q
        q=temp1;
        
        %%%%%%%% LOAD %%%%%%%%
        % changing this to raw_load_cham2 to load the new chameleon files
        dummy = num2str(cast + 10000);
        %    load_file = [path_raw cruise_id '_' dummy(2) '.' dummy(3:5)];
        load_file = fullfile(path_raw,[ cruise_id '_' dummy(2) '.' dummy(3:5)]);
        [data head]=raw_load_cham2(load_file);
        
        %%%%%%%% CALIBRATE %%%%%%%%
        % calibrate the raw voltages into useful data using the calibration
        % coefficients saved in the header OR new calibration coefficients
        % which are now defined in MODIFY_HEADER_**
        cali_eq14;
        
        % make a new modified structure that I will use for chipod method
        cal2=struct();
        
        %~~ check if TP,p,c,t sampled at same rate; if not, interpolate to same
        % length; I want all to be same length as TP
        irTP=head.irep.TP;
        if head.irep.P~=irTP
            clear P2
            P2=nan*ones(length(cal.TP),1);
            P2(1:head.irep.TP:end)=cal.P;
            P2=NANinterp(P2);
            cal2.P=P2;
            cal2.TP=cal.TP;
        end
        
        if head.irep.T1~=irTP
            cal2.T1=interp1(cal.P,cal.T1,cal2.P);
        end
        
        if head.irep.SAL~=irTP
            cal2.SAL=interp1(cal.P,cal.SAL,cal2.P);
        end
        
        if head.irep.FALLSPD~=irTP
            cal2.FALLSPD=interp1(cal.P,cal.FALLSPD,cal2.P);
        end
        %~~
        
        %%%%%%%% SAVE INDIVIDUAL CASTS %%%%%%%%
        temp=num2str(q.script.num+10000);
        fn=[q.script.prefix '_' temp(2:5)];
        
        cal.MakeInfo  = ['Made by AP ' datestr(now) ' w/ ProcessEq14Cham_AP.m'];
        cal2.MakeInfo = ['Made by AP ' datestr(now) ' w/ ProcessEq14Cham_AP.m'];
        save( fullfile( path_save,fn) ,'head','cal','cal2')
        %
        
    end % try
    
end % cast #

%%