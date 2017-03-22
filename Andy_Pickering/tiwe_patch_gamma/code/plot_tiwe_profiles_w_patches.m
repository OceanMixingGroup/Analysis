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

cnum=3102

% patch options
patch_size_min = 1
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
cal.sgth_sm=sw_pden(cal.SAL_sm,cal.T1,cal.P,0);

% load the patches for this *profile*
% if merge_patches==1
%     load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '_merged_minsep_' num2str(min_sep*100) '.mat']) )
% else
%     load(fullfile(save_dir_patch,ot_dir,[project_short '_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_cnum_' num2str(cnum) '.mat']) )
% end

patches = load_tiwe_patches_profile(cnum, patch_size_min,...
    usetemp, merge_patches, min_sep)

%% Plot T,S,sgth profiles

figure(1);clf
set(gcf,'Name',['TIWE Cast ' num2str(cnum)])
agutwocolumn(0.7)
wysiwyg

ax1 = subplot(1,3,1);
plot(cal.T1,cal.P,'linewidth',2)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('T')
ylabel('P [db]','fontsize',16)

ax2 = subplot(1,3,2);
plot(cal.SAL,cal.P,'linewidth',2)
hold on
plot(cal.SAL_sm,cal.P)
%plot(medfilt1(cal.SAL,4),cal.P)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('S','fontsize',16)

ax3 = subplot(1,3,3);
plot(cal.sgth,cal.P,'linewidth',2)
axis ij
grid on
axis tight
ylim([0 200])
xlabel('sgth','fontsize',16)

linkaxes([ax1 ax2 ax3],'y')

%%

tiwe_patches_paths
figdir=fullfile(analysis_dir,project_long,'figures')
print( fullfile( figdir, ['tiwe_cast_' num2str(cnum) '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_T_S_sgth']), '-dpng')

%%
figure(2);clf
agutwocolumn(1)
wysiwyg

rr=1
cc=2

ax1 = subplot(rr,cc,1);
hf=ShadePatchDepths(patches.p1,patches.p2,[14 30])
hold on
plot(cal.T1,cal.P,'k','linewidth',2)
axis ij
axis tight
ylim([0 200])
xlabel('T')
grid on

ax2 = subplot(rr,cc,2);
hf=ShadePatchDepths(patches.p1,patches.p2,[1022 1026])
hold on
plot(cal.sgth,cal.P,'k','linewidth',2)
axis ij
axis tight
ylim([0 200])
xlabel('sgth')
grid on

linkaxes([ax1 ax2],'y')

%%

tiwe_patches_paths
figdir=fullfile(analysis_dir,project_long,'figures')
print( fullfile( figdir, ['tiwe_cast_' num2str(cnum) '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_T_sgth_patches']), '-dpng')

%% Make plots zoomed-in on one patch

ip=10
p1=patches.p1(ip)
p2=patches.p2(ip)

igp=isin(cal.P,[p1 p2]);
igp2=isin(cal.P,[p1-1 p2+1]);


figure(1);clf
set(gcf,'Name',['TIWE Cast ' num2str(cnum)])
agutwocolumn(0.7)
wysiwyg

ax1 = subplot(1,3,1);
plot(cal.T1(igp2),cal.P(igp2),'linewidth',2, 'color',0.5*[1 1 1])
hold on
plot(cal.T1(igp),cal.P(igp),'k','linewidth',2)
axis ij
grid on
axis tight
xlabel('T')
ylabel('P [db]','fontsize',16)
hline(p1,'k--')
hline(p2,'k--')

ax2 = subplot(1,3,2);
plot(cal.SAL(igp2),cal.P(igp2),'linewidth',2, 'color',0.5*[1 1 1])
hold on
plot(cal.SAL(igp),cal.P(igp),'k','linewidth',2)
plot(cal.SAL_sm(igp2),cal.P(igp2),'g--')
plot(cal.SAL_sm(igp),cal.P(igp),'g--')
axis ij
grid on
axis tight
xlabel('S','fontsize',16)
ytloff
hline(p1,'k--')
hline(p2,'k--')

ax3 = subplot(1,3,3);
plot(cal.sgth(igp2),cal.P(igp2),'linewidth',2, 'color',0.5*[1 1 1])
hold on
plot(cal.sgth(igp),cal.P(igp),'k','linewidth',2)
plot(cal.sgth_sm(igp),cal.P(igp),'g--','linewidth',2)
axis ij
grid on
axis tight
xlabel('sgth','fontsize',16)
ytloff
hline(p1,'k--')
hline(p2,'k--')

linkaxes([ax1 ax2 ax3],'y')

%

tiwe_patches_paths
figdir=fullfile(analysis_dir,project_long,'figures')
figname=['tiwe_cast_' num2str(cnum) '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_T_S_sgth_ZOOM_patch_' num2str(ip)]
print( fullfile( figdir, figname), '-dpng')
 
%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
h1=scatter(cal.T1(igp),cal.SAL(igp),'o')
hold on
h2=scatter(cal.T1(igp),cal.SAL_sm(igp),'d')
grid on
xlabel('T')
ylabel('S')
legend([h1 h2],'S','Ssm','location','best')


subplot(222)
h1=scatter(cal.T1(igp),cal.sgth(igp),'o')
hold on
h2=scatter(cal.T1(igp),cal.sgth_sm(igp),'d')
grid on
xlabel('T')
ylabel('sgth')
legend([h1 h2],'no sm','smooth','location','best')


subplot(223)
h1=scatter(cal.SAL(igp),cal.sgth(igp),'o')
hold on
h2=scatter(cal.SAL_sm(igp),cal.sgth_sm(igp),'d')
grid on
xlabel('S')
ylabel('sgth')
legend([h1 h2],'no sm','smooth','location','best')


tiwe_patches_paths
figdir=fullfile(analysis_dir,project_long,'figures')
figname=['tiwe_cast_' num2str(cnum) '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_T_S_sgth_scatter_patch_' num2str(ip)]
print( fullfile( figdir, figname), '-dpng')
%%

figure(2);clf
agutwocolumn(1)
wysiwyg

rr=1
cc=2

ax1 = subplot(rr,cc,1);
%hf=ShadePatchDepths(patches.p1(igp),patches.p,[14 30])
%hold on
plot( cal.T1(igp2), cal.P(igp2), 'color',0.5*[1 1 1],'linewidth',2)
hold on
plot( cal.T1(igp), cal.P(igp), 'k','linewidth',2)
axis ij
axis tight
ylim([p1-1 p2+1])
xlabel('T')
grid on

%
ax2 = subplot(rr,cc,2);
plot( cal.sgth(igp2), cal.P(igp2), 'color',0.5*[1 1 1],'linewidth',2)
hold on
plot(cal.sgth(igp),cal.P(igp),'k','linewidth',2)
axis ij
axis tight
ylim([p1-1 p2+1])
xlabel('sgth')
grid on

linkaxes([ax1 ax2],'y')

%%
