%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotKmaxKbVsEps.m
%
% * Makes plot for chipod methods paper *
%
% Plot the % of theoretical K_b resolved for various values of epsilon
%
%------------------
% 04/25/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

figure(1);clf
agutwocolumn(0.6)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

nu=1e-6;
b_freq=(10.^(-2:.1:3.5))';
fspd=1
tdif=1.5e-7
chi=1e-8
%q=7
fmax=7

eps=[1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4]%
kmax=fmax/fspd
kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
%b_spec= kraichnan(nu,b_freq/fspd,kb,tdif,chi,q)/fspd;

rat1=fmax./kb*100;
rat2=fmax/0.75./kb*100;
rat3=fmax/0.5./kb*100;
rat4=fmax/0.25./kb*100;

h1=semilogx(eps,rat1,'o-','linewidth',2 )
hold on
h2=semilogx(eps,rat2,'o-','linewidth',2 )
h3=semilogx(eps,rat3,'o-','linewidth',2 )
h4=semilogx(eps,rat4,'o-','linewidth',2 )
grid on
ylim([0 200])

pos=get(gca,'position')
set(gca,'position',pos.*[1.2 1.2 1 1])
xlabel('\epsilon [Wkg^{-1}]','fontsize',18)
ylabel('(k_{max}/k_b) *100','fontsize',16)

legend([h1 h2 h3 h4],'1m/s','0.75m/s','0.5m/s','0.25m/s')

%%
SetPaperFigPath
print(fullfile(figdir,'kbratVsEps'),'-dpng')

%%