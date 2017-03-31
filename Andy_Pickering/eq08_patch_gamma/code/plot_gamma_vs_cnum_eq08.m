%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_vs_cnum_eq08.m
%
% Plot gamma estimated for EQ08 patches vs cast #
%
%
%
%-------------
% 3/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.4 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
merge_patches=0
min_sep=0;

eq08_patches_paths

datdir = fullfile( analysis_dir, project_long, 'data')

% load my patches
%load(fullfile(datdir,['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )
patches = load_patches_comb(project_short,patch_size_min,usetemp,merge_patches,min_sep)

depth_range = [60 200]
id = find(patches.p1>depth_range(1) & patches.p2<depth_range(2));

figure(1);clf
agutwocolumn(1)
wysiwyg


%~~~~~~~~~~~~~~~~
subplot(2,2,1)
plot(log10(patches.gam_bin),patches.cnum,'.','color',0.65*[1 1 1])
freqline(log10(0.2))
grid on
xlabel('log_{10}\gamma_{\chi\epsilon} bin','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-2.5 1])
title([project_short ' patches - minOT=' num2str(100*patch_size_min) 'cm, ' ...
    num2str(depth_range(1)) ' - ' num2str(depth_range(2)) 'm'])


% compute median gamma for each cast#
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_bin(id));
end
hold on
plot(log10(gam_md),cnums,'k.','markersize',15)




%~~~~~~~~~~~~~~~~
subplot(2,2,2)

plot(log10(patches.gam_line),patches.cnum,'.','color',0.65*[1 1 1])
freqline(log10(0.2))
grid on
xlabel('log_{10}\gamma_{\chi\epsilon} line','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-2.5 1])


% compute median gamma for each cast#
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_line(id));
end
hold on
plot(log10(gam_md),cnums,'k.','markersize',15)



%~~~~~~~~~~~~~~~~
subplot(2,2,3)
plot(log10(patches.gam_bulk),patches.cnum,'.','color',0.65*[1 1 1])
freqline(log10(0.2))
grid on
xlabel('log_{10}\gamma_{\chi\epsilon} bulk','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-2.5 1])


% compute median gamma for each cast#
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_bulk(id));
end
hold on
plot(log10(gam_md),cnums,'k.','markersize',15)



%~~~~~~~~~~~~~~~~
subplot(2,2,4)
plot(log10(patches.gam_line_fit),patches.cnum,'.','color',0.65*[1 1 1])
freqline(log10(0.2))
grid on
xlabel('log_{10}\gamma_{\chi\epsilon} linefit','fontsize',16)
ylabel('cnum','fontsize',16)
xlim([-2.5 1])


% compute median gamma for each cast#
cnums=1:3200;
gam_md=nan*ones(size(cnums));
for i=1:length(cnums)
    clear id
    id=find(patches.cnum==cnums(i));
    gam_md(i)=nanmedian(patches.gam_line_fit(id));
end
hold on
plot(log10(gam_md),cnums,'k.','markersize',15)



%%

fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)...
    '_gam_vs_cnum_depths_' num2str(depth_range(1)) '_' num2str(depth_range(2))]
print(fullfile(fig_dir,fname),'-dpng')


%%