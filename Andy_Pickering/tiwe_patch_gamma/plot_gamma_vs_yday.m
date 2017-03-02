%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_yday.m
%
% Plot gamma estimated for TIWE patches vs yday
%
%-------------
% 2/24/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.15 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

saveplot=0

% load my patches
load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

figure(1);clf
agutwocolumn(0.75)
wysiwyg

plot(log10(patches.gam_line),patches.yday,'.')
freqline(log10(0.2))
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('yday','fontsize',16)
xlim([-5 3])
title(['TIWE patches - minOT=' num2str(100*patch_size_min) 'cm' ])

ydays=308:329;
gam_md=nan*ones(size(ydays));
for i=1:length(ydays)
    clear id
    id=find(patches.yday==ydays(i));
    gam_md(i)=nanmedian(patches.gam_line(id));
end
hold on
plot(log10(gam_md),ydays,'ro')

if saveplot==1
fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gam_vs_yday']
print(fullfile(fig_dir,fname),'-dpng')
end

%%