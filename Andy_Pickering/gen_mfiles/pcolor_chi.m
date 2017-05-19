function h = pcolor_chi(chipod,cham) 

h=figure(1);clf
agutwocolumn(1)
wysiwyg

rr=4 ;
cc=1 ;

ax1 = subplot(rr,cc,1);
ezpc(cham.cnum,cham.P,log10(cham.chi))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi chameleon')

ax2 = subplot(rr,cc,2);
ezpc(chipod.cnum,chipod.P,log10(chipod.chi))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-11 -4])
colorbar
title('log_{10} \chi \chi-pod')

ax3 = subplot(rr,cc,3);
ezpc(chipod.cnum,chipod.P,real(log10(cham.N2)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-6 -2])
colorbar
ylabel('P [db]')
title('log_{10} N^2')

ax4 = subplot(rr,cc,4);
ezpc(chipod.cnum,chipod.P,real(log10(cham.Tz)))
hold on
%plot(zml_cnum,zml,'k')
%hline(80,'k--')
caxis([-4 -0])
colorbar
ylabel('P [db]')
xlabel('cast #')
title('log_{10} dT/dz')

linkaxes([ax1 ax2 ax3 ax4])

%%