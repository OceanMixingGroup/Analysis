function h = chi_vs_eps_normalized_plot(eps,chi,N2,Tz)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot 2D histogram of chi vs epsilon, normalized such that the slope is
% proportional to gamma
%
%-------------
% A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

h=histogram2(  real(log10(eps./N2)),log10(chi./(Tz.^2)),80,'DisplayStyle','tile');
grid on
hold on
xvec=linspace(1e-7,1e-1,100);
h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
ylim([-8 -1])
xlim([-6.5 -1])
ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','south','orientation','horizontal')

%%