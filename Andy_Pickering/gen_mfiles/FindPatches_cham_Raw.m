function [] = FindPatches_cham_Raw(project_name,patch_size_min,usetemp, cnums_to_do)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FindPatches_cham_Raw.m
%
% * general version for any experiment *
%
% Find patches (overturns) in chameleon profiles, using raw (not
% binned/averaged) data.
%
%
% Dependencies:
%   - compute_overturns_discrete_AP.m
%
%-----------------
% 3/24/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

addpath(fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering',[project_name '_patch_gamma'],'code'))

eval([project_name '_patches_paths'])

save_dir = fullfile( save_dir_patch,['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)],'raw')
ChkMkDir(save_dir)

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/

% loop through each cast
warning off
hb=waitbar(0,['Finding Patches For ' project_name]);

ic=0;
for cnum= cnums_to_do;
    ic=ic+1;
    waitbar(ic/length(cnums_to_do),hb)
    
    try
        
        close all
        clear cal cal2 head patch_data
        
        patch_data=[];
        
        % Load the data for this cast
        %load(fullfile(save_dir_cal,['eq14_' sprintf('%04d',cnum) '.mat']))
        eval( ['cal = load_cal_' project_name '(cnum);'])
        
        %        cal=cal2; clear cal2
        
        clear s t p lat
        %s=smooth( cal.SAL(1:end-1), 20 ); % (end-1) b/c last 2 values are same;
        s = cal.SAL(1:end-1);
        t = cal.T1(1:end-1);
        p = cal.P(1:end-1) ;
        %
        %         clear idot lat1 lat2
        %         idot=strfind(head.lat.start,'.');
        %         lat1=str2num(head.lat.start(1:idot-3));
        %         lat2=str2num(head.lat.start(idot-2:end))/60;
        %         lat=nanmean([lat1 lat2]);
        
        clear Params
        Params.lat=cal.lat;
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
