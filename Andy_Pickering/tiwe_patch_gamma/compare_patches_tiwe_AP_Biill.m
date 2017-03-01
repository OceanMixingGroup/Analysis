%~~~~~~~~~~~~~~~~~~~
%
% compare_patches_tiwe_AP_Bill.m
%
% Compare my patch estimates for TIWE to the results of Bill's analysis. I
% don't have the code he used to produce it, just the resulting patches.
%
%
%------------
% 2/21/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.15 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density
datdir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data'

% load my patches
load(fullfile(datdir,['tiwe_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

% load Bills patches
load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/events_TIWE.mat')

day_range=[307 329]% all profiles
%day_range=[324 327]% ydays used in Smyth etal
id=find(patches.yday>=day_range(1) & patches.yday<=day_range(2));

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
h1=histogram(real(log10(patches.n2_line(id))),'Normalization','pdf','Edgecolor','none');
hold on
h2=histogram(log10(A.N2(:)),h1.BinEdges,'Normalization','pdf','Edgecolor','none')
xlabel('log_{10}[N^2]','fontsize',16)
ylabel('pdf','fontsize',16)
xlim([-6 -2])
grid on
legend([h1 h2],'AP','Bill','location','best')

subplot(222)
hap=histogram(real(log10(patches.dtdz_line(id))),'Normalization','pdf','Edgecolor','none');
hold on
hbill=histogram(real(log10(A.tgrad(:))),hap.BinEdges,'Normalization','pdf','Edgecolor','none');
xlabel('log_{10}[T_z]','fontsize',16)
ylabel('pdf','fontsize',16)
xlim([-4 0])
grid on

subplot(223)
h1=histogram(log10(patches.chi(id)),'Normalization','pdf','Edgecolor','none');
hold on
histogram(real(log10(A.chi(:))),'Normalization','pdf','Edgecolor','none');
xlabel('log_{10}[\chi]')
ylabel('pdf')
xlim([-13 -3])
grid on

subplot(224)
h1=histogram(log10(patches.eps(id)),'Normalization','pdf','Edgecolor','none');
hold on
h2=histogram(real(log10(A.eps(:))),'Normalization','pdf','Edgecolor','none');
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('pdf','fontsize',16)
grid on
legend([h1 h2],'AP','Bill','location','best')
%%

fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_n2_tz_chi_eps_apvsbill_hist_yday_' num2str(day_range(1)) '_' num2str(day_range(2)) ]
print(fullfile(fig_dir,fname),'-dpng')


%% compute gamma from Bill's patches

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles

gam_bill = ComputeGamma(A.N2,A.tgrad,A.chi,A.eps);

figure(2);clf
hbill=histogram(log10(gam_bill(:)),'Normalization','pdf','EdgeColor','none')
hold on
hap=histogram(log10(patches.gam_line(id)),'Normalization','pdf','EdgeColor','none')
freqline(log10(0.2))
xlim([-3 2])
ylim([0 1.2])
grid on
xlabel('log_{10}[\Gamma]')
ylabel('pdf')
legend([hbill hap],'Bill','AP','location','best')
title(['tiwe patches, minOT=' num2str(patch_size_min) 'm, yday ' num2str(day_range(1)) '-' num2str(day_range(2))])

%%
fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
fname=['tiwe_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_gam_apvsbill_hist_yday_' num2str(day_range(1)) '_' num2str(day_range(2))]
print(fullfile(fig_dir,fname),'-dpng')

%%

figure(3);clf
agutwocolumn(0.5)
wysiwyg
h1=histogram(A.Lt(:),'Normalization','pdf');
hold on
h2=histogram(patches.Lt(:),h1.BinEdges,'Normalization','pdf');
legend([h1 h2],'Bill','AP')
grid on
xlabel('L_t')

%%
