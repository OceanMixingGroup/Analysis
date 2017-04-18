function h = PlotEpsProfileCompare_eq14(cnum,Params,dz)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eq14_patches_paths

% regular chi-pod method on binned data
load( fullfile(path_chipod_bin,['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
chb=avg;clear avg

% chamelon data
load(fullfile(path_cham_avg,['EQ14_' sprintf('%04d',cnum) '_avg.mat']))

h = figure(1);clf
set(h,'Name',['eq14 profile ' num2str(cnum)])

ax1 = subplot(121);
hcham = plot(log10(avg.EPSILON),avg.P,'k','linewidth',2);
hold on
hchi = plot(log10(chb.eps1),chb.P,'.','markersize',15,'color',0.5*[1 1 1]);
axis ij
grid on
xlim([-11 -4])
ylim([0 200])
legend([hcham hchi ],'cham','\chi pod','location','best')
xlabel('log_{10}[\epsilon]')
ylabel('P [db]')

[bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
[bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);

% plot 10m binned profiles also
ax2 = subplot(122);
h1 = plot(log10(bin1),z1,'ko-','linewidth',2);
hold on
h2 = plot(log10(bin2),z2,'d-','color',0.5*[1 1 1],'linewidth',2);
axis ij
grid on
xlim([-11 -4])
ylim([0 200])
legend([h1 h2 ],'cham','\chi pod','location','best')
xlabel('log_{10}[\epsilon]')
ylabel('P [db]')
title([num2str(dz) 'm binned'])

linkaxes([ax1 ax2])

return

%%
