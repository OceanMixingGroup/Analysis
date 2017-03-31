function cal = load_cal_tiwe(cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% cal = load_cal_tiwe(cnum)
%
% Load high-res/raw Chameleon cast from tiwe ('cal'); this is calibrated but not
% averaged/binned
%
% These are made in Process_tiwe_rawprofiles.m
%
% 3/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

tiwe_patches_paths
load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )
cal=cal2 ; clear cal2

clear yday
cal.yday=str2num(head.starttime(end-5:end));

% get latitude for profile
clear idot lat1 lat2
idot=strfind(head.lat.start,'.');
lat1=str2num(head.lat.start(1:idot-3));
lat2=str2num(head.lat.start(idot-2:end))/60;
cal.lat=nanmean([lat1 lat2]);

%%