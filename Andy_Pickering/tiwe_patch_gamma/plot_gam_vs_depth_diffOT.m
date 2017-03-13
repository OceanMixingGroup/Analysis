%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_gam_vs_depth_diffOT.m
%
% Plot gamma vs depth for different overturn params
%
%
%-------------
% 3/10/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

day_range=[324 327]% ydays in Smyth et al

p1=load_tiwe_patches_comb(0.4, 1, 0, .15) ;
id1 = find(p1.yday>=day_range(1) & p1.yday<=day_range(2) );

p2=load_tiwe_patches_comb(0.4, 0, 0, .15) ;
id2 = find(p2.yday>=day_range(1) & p2.yday<=day_range(2) );

% Plot gamma vs depth

figure(1);clf
agutwocolumn(0.7)
wysiwyg

subplot(121)
histogram2( log10(p1.gam_line(id1)),p1.p1(id1),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
axis ij
xlim([-3 1.5])
ylim([0 200])
xlabel('log_{10}\gamma','fontsize',16)
ylabel('P','fontsize',16)
title('OT from temp')

subplot(122)
histogram2( log10(p2.gam_line(id2)),p2.p1(id2),50,'DisplayStyle','tile')
freqline(log10(0.2))
grid on
axis ij
xlim([-3 1.5])
ylim([0 200])
xlabel('log_{10}\gamma','fontsize',16)
ylabel('P','fontsize',16)
title('OT from dens')

%%