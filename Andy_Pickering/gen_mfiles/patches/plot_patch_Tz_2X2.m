function h=plot_patch_Tz_2X2(project_name,patch_size_min,usetemp,...
    merge_patches,min_sep,depth_range)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patch_Tz_2X2.m
%
% General version of plotting function for patches/gamma analyses
% 
%
%
%------------------
% 3/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

saveplots=0


eval([project_name '_patches_paths'])

datdir = fullfile(analysis_dir,project_long,'data');


% load my patches
clear patches
patches = load_patches_comb(project_name, patch_size_min, usetemp, merge_patches, min_sep)

%%

id = find( patches.p1>depth_range(1) & patches.p2<depth_range(2) );
id2 = find( patches.p1>depth_range(1) & patches.p2<depth_range(2) & patches.R2>0.5 );
length(id)

%%

xl=[-4 0]

h=figure;clf
agutwocolumn(1)
wysiwyg
set(gcf,'name',['min OT ' num2str(patch_size_min) 'm'])

ax1 = subplot(221) ;
h1=histogram(real(log10(patches.dtdz_bin(id))),'Normalization','pdf', 'EdgeColor','none');
hold on
histogram(real(log10(patches.dtdz_bin(id2))),'Normalization','pdf', 'EdgeColor','none');
xlim(xl)
freqline(-2);
grid on
xlabel('log_{10}[T_z bin]','fontsize',16)
ylabel('pdf','fontsize',16)
      
ax2 = subplot(222) ;
h2=histogram(real(log10(patches.dtdz_line(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
hold on
histogram(real(log10(patches.dtdz_line(id2))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim(xl)
freqline(-2);
grid on
xlabel('log_{10}[T_z line]','fontsize',16)
ylabel('pdf','fontsize',16)

ax3 = subplot(223) ;
h3=histogram(real(log10(patches.dtdz_bulk(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
hold on
histogram(real(log10(patches.dtdz_bulk(id2))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim(xl)
freqline(-2);
grid on
xlabel('log_{10}[T_z bulk]','fontsize',16)
ylabel('pdf','fontsize',16)

ax4 = subplot(224) ;
% h4=histogram(real(log10(patches.dtdz_line_fit(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
% hold on
% histogram(real(log10(patches.dtdz_line_fit(id2))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
% xlim([-6.5 -2])
% freqline(log10(0.2));
% grid on
% xlabel('log_{10}[T_z line-fit]','fontsize',16)
% ylabel('pdf','fontsize',16)
%legend([h1 h2 h3 h4],'bin','line','bulk','line-fit','location','best')

linkaxes([ax1 ax2 ax3 ax4])


%%

if saveplots==1
    
    if merge_patches==1
        fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist2X2_merged_minsep_' num2str(min_sep*100) '_depth_' num2str(depth_range(1)) '_' num2str(depth_range(2)) 'm' ]
    else
        fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist2X2_depth_' num2str(depth_range(1)) '_' num2str(depth_range(2)) 'm' ]
    end
    
%    print(fullfile(fig_dir,fname),'-dpng')
end

%%