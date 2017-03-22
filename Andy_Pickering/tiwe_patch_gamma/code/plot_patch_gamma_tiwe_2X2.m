function h=plot_patch_gamma_tiwe_2X2(patch_size_min,usetemp,...
    merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patch_gamma_tiwe.m
%
%
%------------------
% 2/20/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

saveplots=1

tiwe_patches_paths

datdir=fullfile(analysis_dir,project_long,'data') ; 

% load my patches
clear patches
if merge_patches==1
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma_merged_minsep_' num2str(min_sep*100)  '.mat']) )
else
    load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']) )
end

%day_range=[307 329]% all profiles
day_range=[324 327]% ydays in Smyth et al
depth_range= [60 200]
id = find(patches.yday>=day_range(1) & patches.yday<=day_range(2) & patches.p1>depth_range(1) & patches.p2<depth_range(2) );

length(id)
%

h=figure;clf
agutwocolumn(1)
wysiwyg
set(gcf,'name',['min OT ' num2str(patch_size_min) 'm'])

ax1 = subplot(221) ;
h1=histogram(real(log10(patches.gam_bin(id))),'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
freqline(log10(0.2));
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon} bin]','fontsize',16)
ylabel('pdf','fontsize',16)
%title(['min OT ' num2str(patch_size_min) 'm'])
      
ax2 = subplot(222) ;
h2=histogram(real(log10(patches.gam_line(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
freqline(log10(0.2));
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon} line]','fontsize',16)
ylabel('pdf','fontsize',16)

ax3 = subplot(223) ;
h3=histogram(real(log10(patches.gam_bulk(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
freqline(log10(0.2));
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon} bulk]','fontsize',16)
ylabel('pdf','fontsize',16)

ax4 = subplot(224) ;
h4=histogram(real(log10(patches.gam_line_fit(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
freqline(log10(0.2));
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon} line-fit]','fontsize',16)
ylabel('pdf','fontsize',16)
%legend([h1 h2 h3 h4],'bin','line','bulk','line-fit','location','best')

linkaxes([ax1 ax2 ax3 ax4])
if merge_patches==1
%    title(['min OT ' num2str(patch_size_min) 'm, yday ' num2str(day_range(1)) '-' num2str(day_range(2)) ', ' num2str(depth_range(1)) ':' num2str(depth_range(2)) ' db, merged, N=' num2str(length(id))])
else
 %   title(['min OT ' num2str(patch_size_min) 'm, yday ' num2str(day_range(1)) '-' num2str(day_range(2)) ', ' num2str(depth_range(1)) ':' num2str(depth_range(2)) ' db, N=' num2str(length(id))])
end

hold on
% plot(nanmedian(log10(patches.gam_bin(id))),0,'d')
% plot(nanmedian(log10(patches.gam_line(id))),0,'d')
% plot(nanmedian(log10(patches.gam_line_fit(id))),0,'d')
% plot(nanmedian(log10(patches.gam_bulk(id))),0,'d')
clear patches
%%

if saveplots==1
    
    if merge_patches==1
        fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist2X2_yday_' num2str(day_range(1)) '_' num2str(day_range(2)) '_merged_minsep_' num2str(min_sep*100)]
    else
        fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist2X2_yday_' num2str(day_range(1)) '_' num2str(day_range(2)) ]
    end
    
    print(fullfile(fig_dir,fname),'-dpng')
end
%%