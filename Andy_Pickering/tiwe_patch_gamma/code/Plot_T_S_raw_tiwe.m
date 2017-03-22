%~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_T_S_raw_tiwe.m
%
%
% 3/8/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

tiwe_patches_paths

T=[];
S=[];
P=[];

for cnum=2836:50:3711
    try
        clear cal cal2 head
        load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )
        
        T = [T(:) ; cal2.T1(:)];
        S = [S(:) ; cal2.SAL(:)];
        P = [P(:) ; cal2.P(:) ];
    end
end

%%

figure(1);clf
agutwocolumn(0.8)
wysiwyg

ax1=subplot(211);
scatter(T,S,'.','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
grid on
ylabel('S')

% also plot for just 60-200m
idz=isin(P,[60 200]);

ax2=subplot(212);
h=scatter(T(idz),S(idz),'.','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
grid on
xlabel('T')
ylabel('S')

linkaxes([ax1 ax2])

%%

figdir=fullfile(analysis_dir,project_long,'figures')
print(fullfile(figdir,'tiwe_T_S_scatter'),'-dpng')

%%


clear ; close all

tiwe_patches_paths

cnum=3000
clear cal cal2 head
load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )

figure(2);clf
scatter(cal2.T1,cal2.SAL,'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
idz=isin(cal2.P,[60 200]);
hold on
scatter(cal2.T1(idz),cal2.SAL(idz),'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1)
grid on

%%

clear ; close all

tiwe_patches_paths

cnum=3102
clear cal cal2 head
load( fullfile( save_dir_cal, ['tw91' sprintf('%04d',cnum) '_raw.mat'] ) )

figure(1);clf

ax1 = subplot(131)
plot(cal2.T1,cal2.P,'linewidth',1);
axis ij
axis tight
ylim([0 200])
grid on
xlabel('T')

ax2 =subplot(132);
plot(cal2.SAL,cal2.P,'linewidth',1);
axis ij
axis tight
ylim([0 200])
xlabel('S')
grid on

addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/
cal2.sgth=sw_pden(cal2.SAL,cal2.T1,cal2.P,0);

ax3 = subplot(133);
plot(cal2.sgth,cal2.P,'linewidth',1)
axis ij
axis tight
ylim([0 200])
xlabel('sgth')
grid on

linkaxes([ax1 ax2 ax3],'y')

ylim([147 160])

%%