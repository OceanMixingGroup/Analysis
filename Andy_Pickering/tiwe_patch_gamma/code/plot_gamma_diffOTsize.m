%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gamma_diffOTsize.m
%
% Plot gamma from tiwe patches for different min OT sizes
%
%---------
% 3/10/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~
%% compare gamma for different OT size

clear ; close all

usetemp=0

patch_size_min1 = 0.4
p1=load_tiwe_patches_comb(patch_size_min1, usetemp, 0, .15) ;
id1 = find(p1.p1>60 & p1.p2<200);

patch_size_min2 = 1
p2=load_tiwe_patches_comb(patch_size_min2, usetemp, 0, .15) ;
id2 = find(p2.p1>60 & p2.p2<200);


figure(1);clf
agutwocolumn(0.8)
wysiwyg

ax1 = subplot(211) ;
h1 = histogram( real(log10(p1.gam_bin(id1))),'Normalization','pdf','Edgecolor','none');
hold on
h2 = histogram( real(log10(p2.gam_bin(id2))),'Normalization','pdf','Edgecolor','none');
xlim([-3 2])
freqline(log10(0.2))
grid on
xlabel('log_{10}[\gamma]','fontsize',16)
ylabel('pdf')
title(['  \gamma bin'])
legend([h1 h2],num2str(patch_size_min1), num2str(patch_size_min2))

ax2 = subplot(212);
h1 = histogram( real(log10(p1.gam_line(id1))),'Normalization','pdf','Edgecolor','none');
hold on
h2 = histogram( real(log10(p2.gam_line(id2))),'Normalization','pdf','Edgecolor','none');
xlim([-4 2])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]','fontsize',16)
ylabel('pdf')
grid on
title([' \gamma line'])
legend([h1 h2],num2str(patch_size_min1), num2str(patch_size_min2))

%%

tiwe_patches_paths
figname=['tiwe_usetemp_' num2str(usetemp) '_gam_diffMinOT']
print(fullfile(fig_dir,figname),'-dpng')


%%