function [] = FindPatches_EQ14_Raw(patch_size_min,usetemp, cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FindPatches_EQ14_Raw.m
%
% Find patches (overturns) in EQ14 chameleon profiles, using raw (not
% binned/averaged) data.
%
%
% The raw mat files for each chameleon cast are made w/
% ProcessEq14Cham_AP.m, which was modified from Sally's code so I could
% make files to apply chipod method to. See also ComputeChi_Chameleon_Eq14.m
% They are in the folder /ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/
%
% Dependencies:
%   - compute_overturns_discrete_AP.m
%
%-----------------
% 10/27/16 - A.Pickering - andypicke@gmail.com
% 11/01/16 - AP - Use specific cast #s so we can match up with other data
% later more easily
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%clear ; close all

% patch options
%patch_size_min = 0.15 ; % min patch size
%usetemp   = 1 ;         % 1=use pot. temp, 0= use density

eq14_patches_paths

save_dir = fullfile( save_dir_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)],'raw')
ChkMkDir(save_dir)

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/

% loop through each cast
warning off
hb=waitbar(0,'FindPatches EQ14 Raw');

% only do profiles that are done in chameleon processing
%cnums_to_do=[4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
%    1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];

ic=0;
for cnum= cnums_to_do;
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    
    try
        
        close all
        clear cal cal2 head patch_data
        
        patch_data=[];
        
        % Load the data for this cast
        load(fullfile(save_dir_cal,['eq14_' sprintf('%04d',cnum) '.mat']))
        
        cal=cal2; clear cal2
        
        clear s t p lat
        %s=smooth( cal.SAL(1:end-1), 20 ); % (end-1) b/c last 2 values are same;
        s = cal.SAL(1:end-1);
        t = cal.T1(1:end-1);
        p = cal.P(1:end-1) ;
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        lat=nanmean([lat1 lat2]);
        
        clear Params
        Params.lat=lat;
        Params.plotit=0;
        Params.sigma=1e-5;
        Params.runlmin=0;
        Params.minotsize=patch_size_min;
        Params.usetemp=usetemp;
        
        clear OT
        OT=compute_overturns_discrete_AP(p,t,s,Params);
        
        clear pstarts pstops
        pstarts=OT.pstarts_each;
        pstops=OT.pstops_each;
        
        for i=1:length(pstarts)
            % don't keep patches shallower than 10m depth
            if pstarts(i)>10
                patch_data=[patch_data ; cnum pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) OT.Otnsq_each(i) OT.Lt_each(i) ];
            end
        end
        
        col_names={'cnum','pstart','pstop','dp','otnsq','otLt'};
        
        % save patch data for this profile
        save(fullfile(save_dir,[project_short '_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']),'patch_data','col_names')
        
    end % try
    
end % cnum

delete(hb)
warning on


%%