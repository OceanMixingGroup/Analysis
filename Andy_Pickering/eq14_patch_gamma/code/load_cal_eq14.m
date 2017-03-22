function cal = load_cal_eq14(cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% cal = load_cal_eq14(cnum)
%
% Load high-res/raw Chameleon cast from eq14 ('cal'); this is calibrated but not
% averaged/binned
%
% These are made in xxxx.m
%
%-----------------
% 3/17/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

eq14_patches_paths
%load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )
%cal=cal2 ; clear cal2
load(fullfile(save_dir_cal,['eq14_' sprintf('%04d',cnum) '.mat']))       
cal=cal2 ; clear cal2

% % get latitude for profile
% clear idot lat1 lat2
% idot=strfind(head.lat.start,'.');
% lat1=str2num(head.lat.start(1:idot-3));
% lat2=str2num(head.lat.start(idot-2:end))/60;
% cal.lat=nanmean([lat1 lat2]);

%%