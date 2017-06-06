function h = scatter_chi_eps_chipod_cham(chipod,cham)

%%

h = figure;clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(211) ;
histogram2( log10(cham.chi(:)), log10(chipod.chi(:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-12 -4])
ylim([-12 -4])
xlabel('\chi','fontsize',16)
ylabel('\chi_{\chi}','fontsize',16)

ax2 = subplot(212);
histogram2( log10(cham.eps(:)), log10(chipod.eps(:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-8.5 -4])
ylim([-8.5 -4])
xlabel('\epsilon ','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)

%%