%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patch_gamma_tiwe.m
%
%
%------------------
% 2/20/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1

load( fullfile( '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/',...
    ['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))

%%
day_range=[307 329]% all profiles
day_range=[324 327]% ydays in Smyth et al
id=find(patches.yday>=day_range(1) & patches.yday<=day_range(2));

%

figure(1);clf
agutwocolumn(0.6)
wysiwyg
h1=histogram(real(log10(patches.gam_bin(id))),'Normalization','pdf', 'EdgeColor','none');
hold on
h2=histogram(real(log10(patches.gam_line(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
h3=histogram(real(log10(patches.gam_bulk(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
%ylim([0 1.2])
freqline(log10(0.2))
grid on
xlabel('log_{10}[\Gamma]','fontsize',16)
ylabel('pdf','fontsize',16)
legend([h1 h2 h3],'bin','line','bulk','location','best')

%%

fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist_yday_' num2str(day_range(1)) '_' num2str(day_range(2)) ]
print(fullfile(fig_dir,fname),'-dpng')

%%