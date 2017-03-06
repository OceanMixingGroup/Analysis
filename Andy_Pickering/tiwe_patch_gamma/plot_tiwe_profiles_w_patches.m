%~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_tiwe_profiles_w_patches.m
%
% Script to plot TIWE profiles and patches.
%
%
%----------------
% 3/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cnum=3120

% patch options
patch_size_min = 0.75
usetemp=0
merge_patches=0
min_sep=0.15

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

tiwe_patches_paths
ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)];

addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

% load raw (un-averaged) profile
% Load the data for this cast
load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )
cal=cal2 ; clear cal2

cal.sgth=sw_pden(cal.SAL,cal.T1,cal.P,0);


% load the patches for this profile
if merge_patches==1
    load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
else
    load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
end

%% Plot T,S,sgth profiles

figure(1);clf
agutwocolumn(0.7)
wysiwyg

ax1 = subplot(1,3,1);
plot(cal.T1,cal.P)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('T')

ax2 = subplot(1,3,2);
plot(cal.SAL,cal.P)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('S')

ax3 = subplot(1,3,3);
plot(cal.sgth,cal.P)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('sgth')

linkaxes([ax1 ax2 ax3],'y')

%%
figure(1);clf
agutwocolumn(1)
wysiwyg

rr=1
cc=2

ax1 = subplot(rr,cc,1);
hf=ShadePatchDepths(patches.p1,patches.p2,[14 30])
hold on
plot(cal.T1,cal.P,'k','linewidth',2)
axis ij
ylim([0 200])
xlabel('T')

ax2 = subplot(rr,cc,2);
hf=ShadePatchDepths(patches.p1,patches.p2,[1022 1026])
hold on
plot(cal.sgth,cal.P,'k','linewidth',2)
axis ij
ylim([0 200])
xlabel('sgth')

linkaxes([ax1 ax2],'y')

%%
