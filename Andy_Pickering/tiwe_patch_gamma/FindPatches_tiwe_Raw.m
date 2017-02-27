%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FindPatches_tiwe_Raw.m
%
% Modifed from FindPatches_EQ14_Raw.m
%
% Find patches (overturns) in tiwe chameleon profiles, using raw (not
% binned/averaged) data.
%
%
% Dependencies:
%   - compute_overturns_discrete_AP.m
%
%-----------------
% 2/15/17 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/

tiwe_patches_paths
%
datdir = save_dir_cal 
save_dir = save_dir_patch
ChkMkDir(save_dir)

% patch options
save_data = 1 ;         % save data at the end
patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density


% loop through each cast
warning off
hb=waitbar(0,'working on profiles');

Flist = dir(fullfile(datdir,'*raw.mat'))
%%

for ic= 1:length(Flist)
    
    waitbar(ic/length(Flist),hb)
    clear patch_data
    patch_data=[];

    try
        
        close all
        clear cal cal2 head fname cnum
        
        fname=Flist(ic).name ;
        cnum=str2num(fname(5:8));
                
        % Load the data for this cast
        load(fullfile(datdir,Flist(ic).name))
        
        clear yday
        yday=str2num(head.starttime(end-5:end));
        
        cal=cal2; clear cal2
        
        clear s t p lat
        %s=smooth( cal.SAL(1:end-1), 20 ); % (end-1) b/c last 2 values are same;
        s=cal.SAL(1:end-1);
        t=cal.T1(1:end-1);
        p=cal.P(1:end-1) ;
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        lat=nanmean([lat1 lat2]);
        if isnan(lat)
            lat=0.3;
        end
        
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
                patch_data=[patch_data ; cnum pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) OT.Otnsq_each(i) OT.Lt_each(i) yday ];
            end
        end
                
    end % try
    
    col_names={'cnum','pstart','pstop','dp','otnsq','otLt','yday'};
    
    % save patch data for this profile
    save(fullfile(save_dir,['tiwe_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_cnum_' num2str(cnum) '.mat']),'patch_data','col_names')
    
end % cnum

delete(hb)
warning on

if save_data==1
    savedir = fullfile(analysis_dir, project, 'data' )    
    fname   = ['tiwe_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '.mat'] 
    save( fullfile( savedir,fname), 'patch_data')
end

%%
