%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patch_gamma_eq14.m
%
%
%------------------
% 2/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.4  % min patch size
usetemp = 1

load( fullfile( '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))

cast_range=[1 3100] % all casts
z_range = [0 200]
%z_range = [60 200]
%id=find(patches.cnum>=cast_range(1) & patches.cnum<=cast_range(2));
id=find(patches.cnum>=cast_range(1) & patches.cnum<=cast_range(2) & patches.p1>z_range(1) & patches.p2<z_range(2));

figure(1);clf
agutwocolumn(0.75)
wysiwyg
h1=histogram(real(log10(patches.gam_bin(id))),'Normalization','pdf', 'EdgeColor','none');
hold on
h2=histogram(real(log10(patches.gam_line(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
h3=histogram(real(log10(patches.gam_bulk(id))),h1.BinEdges,'Normalization','pdf', 'EdgeColor','none');
xlim([-3.5 1.5])
freqline(log10(0.2))
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)
legend([h1 h2 h3],'bin','line','bulk','location','best')

title(['EQ14 patches, minOT=' num2str(patch_size_min) ' m, ' num2str(z_range(1)) '-' num2str(z_range(2)) 'db'])

%%

eq14_patches_paths
fname=['eq14_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist_cast_' num2str(cast_range(1)) '_' num2str(cast_range(2)) '_zrange_' num2str(z_range(1)) '_' num2str(z_range(2)) ]
print(fullfile(figdir,fname),'-dpng')

%%