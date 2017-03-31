function h=plot_gamma_vs_epsilon2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
% 3/22/17 - A.Pickering
% 3/29/17 - AP - add depth_range
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

eval([project_name '_patches_paths'])

datdir = fullfile(analysis_dir,project_long,'data');

% load my patches
clear patches
patches = load_patches_comb(project_name, patch_size_min, usetemp, merge_patches, min_sep)

id = find( patches.p1>depth_range(1) & patches.p2<depth_range(2) ) ;

%id=1:length(patches.p1);

%% Plot gamma vs epsilon

h=figure;clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(221) ;
histogram2( real(log10(patches.gam_bin(id))),log10(patches.eps(id)),50,'DisplayStyle','tile');
freqline(log10(0.2))
grid on
xlim([-3.5 1.5])
xlabel('log_{10}[\gamma bin]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)

ax2 = subplot(222) ;
histogram2( log10(patches.gam_line(id)),log10(patches.eps(id)),50,'DisplayStyle','tile');
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma line]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)

ax3 = subplot(223) ;
histogram2( log10(patches.gam_bulk(id)),log10(patches.eps(id)),50,'DisplayStyle','tile');
freqline(log10(0.2))
grid on
xlim([-3 1.5])
xlabel('log_{10}[\gamma bulk]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)

ax4 = subplot(224) ;
histogram2( real(log10(patches.gam_line_fit(id))),log10(patches.eps(id)),50,'DisplayStyle','tile');
freqline(log10(0.2))
grid on
xlim([-3.5 1.5])
xlabel('log_{10}[\gamma line fit]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)

linkaxes([ax1 ax2 ax3 ax4])

if merge_patches==1
    fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_vs_eps2X2_' num2str(depth_range(1)) '_' num2str(depth_range(2)) 'm merged_minsep_' num2str(min_sep*100)]
else
    fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_vs_eps2X2_' num2str(depth_range(1)) '_' num2str(depth_range(2)) 'm' ]
end


print(fullfile(fig_dir,fname),'-dpng')

%%