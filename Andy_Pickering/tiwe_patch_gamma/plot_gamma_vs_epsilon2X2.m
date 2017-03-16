function h=plot_gamma_vs_epsilon2X2(patch_size_min,usetemp,...
    merge_patches,min_sep)
%%

tiwe_patches_paths

datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

% load my patches
clear patches
if merge_patches==1
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) )
else
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) )
end

%day_range=[307 329]% all profiles
day_range=[324 327]% ydays in Smyth et al
%depth_range= [60 200]
%id = find(patches.yday>=day_range(1) & patches.yday<=day_range(2) & patches.p1>depth_range(1) & patches.p2<depth_range(2) );
id = find(patches.yday>=day_range(1) & patches.yday<=day_range(2) );

%% Plot gamma vs epsilon

h=figure;clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(221) ;
histogram2( real(log10(patches.gam_bin)),log10(patches.eps),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma bin]','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

ax2 = subplot(222) ;
histogram2( log10(patches.gam_line),log10(patches.eps),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma line]','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

ax3 = subplot(223) ;
histogram2( log10(patches.gam_bulk),log10(patches.eps),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma bulk]','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

ax4 = subplot(224) ;
histogram2( real(log10(patches.gam_line_fit)),log10(patches.eps),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma line fit]','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

if merge_patches==1
%title(['TIWE patches - minOT=' num2str(100*patch_size_min) 'cm, merged' ])    
else
%title(['TIWE patches - minOT=' num2str(100*patch_size_min) 'cm' ])
end

if merge_patches==1
    fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_vs_eps2X2_' num2str(day_range(1)) '_' num2str(day_range(2)) '_merged_minsep_' num2str(min_sep*100)]
else
    fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_vs_eps2X2_' num2str(day_range(1)) '_' num2str(day_range(2)) ]
end

%fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures';

print(fullfile(fig_dir,fname),'-dpng')

%%