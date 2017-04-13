function h = PlotEpsProfileCompare_eq08(cnum,Params,patch_size_min,...
    usetemp,minR2,dz)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eq08_patches_paths

dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches') ;

% patch N^2,dTdz w/ constant gamma
%load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
%ch=avg;clear avg

% regular chi-pod method on binned data
load(fullfile(path_chipod_bin,['zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['eq08_' sprintf('%04d',cnum) '_avg.mat']) )
chb=avg;clear avg

% chamelon data
%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])
%load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum_filtered.mat')
load(fullfile(path_cham_avg,[project_short '_' sprintf('%04d',cnum) '_avg.mat']) )

%avg=cham;

h = figure(1);clf
agutwocolumn(0.75)
wysiwyg

set(h,'Name',[project_short ' profile ' num2str(cnum)])

ax1 = subplot(121);
h1 = plot(log10(avg.EPSILON),avg.P,'k','linewidth',2);
hold on
h2 = plot(log10(chb.eps1),chb.P,'.','markersize',15,'color',0.5*[1 1 1]);
%h3 = plot(log10(ch.eps1),ch.P,'rp','linewidth',2,'markersize',12);
axis ij
grid on
xlim([-11 -4])
ylim([0 200])
legend([h1 h2 ],'\epsilon cham','\epsilon_{\chi}bin','location','best')
xlabel('log_{10}[\epsilon]')
ylabel('P [db]')


[bin1 z1 Nobs] = binprofile(avg.EPSILON, avg.P, 0, dz, 200,1);
[bin2 z2 Nobs] = binprofile(chb.eps1   , chb.P, 0, dz, 200,1);
%[bin3 z3 Nobs] = binprofile(ch.eps1    , ch.P , 0, dz, 200,0);

% plot 10m binned profiles also
ax2 = subplot(122);
h1 = plot(log10(bin1),z1,'kd-','linewidth',2);
hold on
h2 = plot(log10(bin2),z2,'o-','color',0.5*[1 1 1],'linewidth',2);
%h3 = plot(log10(bin3),z3,'rp-','linewidth',2,'markersize',12);
axis ij
grid on
xlim([-11 -3])
ylim([0 200])
legend([h1 h2 ],'\epsilon cham','\epsilon_{\chi} bin','location','best')
xlabel('log_{10}[\epsilon]')
ylabel('P [db]')

linkaxes([ax1 ax2])

return

%%
