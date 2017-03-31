function h=plot_gamma_vs_yday(patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_yday.m
%
% Plot gamma estimated for TIWE patches vs yday
%
%-------------
% 2/24/17 - A.Pickering
% 3/29/17 - AP - Modify to specify depth range also
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

saveplot=1

% load my patches
clear patches
if merge_patches==1
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) )
else
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) )
end

id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2) );

h=figure;clf
agutwocolumn(0.75)
wysiwyg

ax2=subplot(1,3,[2 3]);
plot(log10(patches.gam_line(id)),patches.yday(id),'.')
freqline(log10(0.2));
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('yday','fontsize',16)
xlim([-3 2.5])
title(['TIWE patches - minOT=' num2str(100*patch_size_min)...
    'cm, depth range ' num2str(depth_range(1)) ' - ' num2str(depth_range(2)) ])

ydays=308:329;
gam_z = patches.gam_line(id);
yday_z=patches.yday(id);
gam_md=nan*ones(size(ydays));
N=nan*ones(size(ydays));
for i=1:length(ydays)
    clear id2
    id2=find(yday_z==ydays(i));
    gam_md(i) = nanmedian(gam_z(id2));
    N(i) = length(id2) ;
end
hold on
plot(log10(gam_md),ydays,'ro');

ax1=subplot(1,3,1);
plot(N,ydays,'o-')
grid on
xlabel('N patches')

linkaxes([ax1 ax2],'y')

if saveplot==1
    fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
    if merge_patches==1
        fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_depths_' num2str(depth_range(1)) '_' num2str(depth_range(2)) '_gam_vs_yday_merged']
    else
        fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_depths_' num2str(depth_range(1)) '_' num2str(depth_range(2)) '_gam_vs_yday']
    end
    print(fullfile(fig_dir,fname),'-dpng')
end

%%