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

load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/TIWE/data/',...
    ['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))

%%

figure(1);clf
agutwocolumn(0.6)
wysiwyg
h1=histogram(real(log10(patches.gam_bin(:))),'Normalization','pdf');
hold on
h2=histogram(real(log10(patches.gam_line(:))),h1.BinEdges,'Normalization','pdf');
h3=histogram(real(log10(patches.gam_bulk(:))),h1.BinEdges,'Normalization','pdf');
xlim([-4 1])
ylim([0 1.2])
freqline(log10(0.2))
grid on
xlabel('log_{10}[\Gamma]','fontsize',16)
ylabel('pdf','fontsize',16)
legend([h1 h2 h3],'bin','line','bulk','location','best')

%%

fig_dir='/Users/Andy/Cruises_Research/ChiPod/TIWE/figures'
fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gammas_hist']
print(fullfile(fig_dir,fname),'-dpng')

%%