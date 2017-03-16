function h=plot_gamma_binned
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%


clear ; close all

% load binned (1m avg) chameleon profiles
load( fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/','tiwe_1mavg_combined.mat') )

% Compute gamma from these values
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
gam=ComputeGamma(cham.N2(:),cham.DTDZ(:),cham.CHI(:),cham.EPSILON(:));

%day_range = [307 329]% all profiles
day_range = [324 327]% ydays used in Smyth etal
depth_range= [60 200]
%id = find(cham.yday>=day_range(1) & cham.yday<=day_range(2) & cham.P>depth_range(1) & cham.P<depth_range(2) );

ib=find(cham.P<60 & cham.P<200);
gam(ib)=nan;

id = find(cham.yday>=day_range(1) & cham.yday<=day_range(2)  ) ;

h = figure ; clf
histogram(real(log10(gam(id))),35,'EdgeColor','none','Normalization','pdf')
freqline(log10(0.2))
%freqline(nanmean(log10(gam(id))),'b--')
xlim([-4 2])
grid on
xlabel('log_{10}[\gamma_{\chi\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)
title(['tiwe 1m bin, yday ' num2str(day_range(1)) '-' num2str(day_range(2))])

%%
fig_dir='/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/figures'
fig_name=['tiwe_avgCombine_gamma_yday_' num2str(day_range(1)) '_' num2str(day_range(2))]
print(fullfile(fig_dir,fig_name),'-dpng')

%%