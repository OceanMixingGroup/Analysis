%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotSpecShapes.m
%
% Plot Nasmyth and Kraichnan spectra for paper
%
% Dependencies:
% - nasmyth.m
% - unv_spec.m
% - kraichnan.m
%
%--------------
% 07/29/16 - A.Pickering - apickering@coas.oregonstate.edu
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=0

figure(1);clf
agutwocolumn(1)
wysiwyg

% Plot Nasmyth spectra for shear
ax1=subplot(211)

freq=(10.^(-2:.1:3.5))';
fspd=1 ;
epsilon=1e-8;
nu=1e-6;
ks=(epsilon/(nu^3))^.25;

[f2,kks]=nasmyth(1000,20);
[unfreq,unspec]=unv_spec(epsilon,nu,kks,f2,fspd);
h1=loglog(unfreq/fspd,unspec,'linewidth',2)

hold on

% Plot again for different epsilon value
clear epsilon ks fs kks
epsilon=1e-5;
ks=(epsilon/(nu^3))^.25;
[f2,kks]=nasmyth(1000,20);
[unfreq,unspec]=unv_spec(epsilon,nu,kks,f2,fspd);
h2=loglog(unfreq/fspd,unspec,'linewidth',2)

% Add legend
hl=legend([h1 h2],'\epsilon=10^{-8}','\epsilon=10^{-5}');
hl.FontSize=16;

xlim([1e-2 1e5])
ylim([1e-7 1e-1])
grid on
xlabel('wavenumber [cpm]','fontsize',18)
ylabel('\Phi_{u_z} [s^{-1}cpm^{-1}] ','fontsize',18)
title('Shear Spectra')
set(gca,'FontSize',15)

text(10^(-1),10^(-3),'\Phi \propto \epsilon','fontsize',16,'fontweight','bold');
annotation('arrow',[0.2 0.2],[0.8 0.9]);
annotation('arrow',[0.2 0.3],[0.8 0.8]);


%~~~~~~~~~

% In 2nd panel, plot Kraichnan spectra for dT/dt

ax2=subplot(212)

b_freq=(10.^(-2:.1:3.5))';
fspd=1;
nu=1e-6;
tdif=1.5e-7;
chi=1e-8;
q=7;
hh=[];

cols=['b' 'r' 'g' 'k'];
ic=1;

% Plot for 2 epsilon values
for eps=[1e-10  1e-5]
    clear kb b_spec h
    kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
    b_spec= kraichnan(nu,b_freq/fspd,kb,tdif,chi,q)/fspd;
    h=loglog(b_freq/fspd,b_spec,cols(ic),'linewidth',2);
    hh=[hh h];
    hold on
    % Also mark batchelor wavenumber
    loglog(kb,1e-7,'d','color',cols(ic),'linewidth',3,'markersize',10)
    ic=ic+1;
end

% Repeat, for a different value of \chi
chi=1e-6;
ic=1;
for eps=[1e-10    1e-5]
    clear kb b_spec h
    kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
    b_spec= kraichnan(nu,b_freq/fspd,kb,tdif,chi,q)/fspd;
    h=loglog(b_freq/fspd,b_spec,'--','color',cols(ic),'linewidth',2);
    loglog(kb,1e-7,'d','color',cols(ic),'linewidth',3,'markersize',10)
    hh=[hh h];
    hold on
    ic=ic+1;
end

% Add legend
hl=legend([hh ],'\epsilon=10^{-10},\chi=10^{-8}','\epsilon=10^{-5},\chi=10^{-8}',...
    '\epsilon=10^{-10},\chi=10^{-6}','\epsilon=10^{-5},\chi=10^{-6}');
hl.FontSize=14;

xlim([1e-2 1e5])
ylim([1e-7 1e-1])
xlabel('wavenumber [cpm]','fontsize',18)
ylabel('\Phi_{T_z} K^2s^{-2}/s^{-1}','fontsize',18)
grid on
set(gca,'FontSize',15)
title('Temperature Gradient Spectra')

% Add some annotations giving scaling of spectra w/ epsilon/chi
text(10^(-1.5),10^(-3),'\Phi \propto \chi,\epsilon^{-1/2}','fontsize',15,'fontweight','bold')
text(10^(3.5),10^(-6),' k_b \propto \epsilon^{1/4}','fontsize',15,'fontweight','bold')
annotation('arrow',[0.8 0.9],[0.1 0.1])

% Mark protion of spectrum we use for chi-pod fits
hh=line([2 2],[1e-10 1e-1]);
hh.LineStyle='--';
hh=line([7 7],[1e-10 1e-1]);
hh.LineStyle='--';

% Plot the k^1/3 slope
xx=linspace(1e-2,1e2,100);
loglog(xx,xx.^(1/3)/1e4,'k--')
linkaxes([ax1 ax2],'x')

% komogorov wave #
%kk= (eps/nu.^3)^(1/4)
%loglog(kk,1e-6,'d','color',cols(ic),'linewidth',3,'markersize',10)

%%
if saveplot==1
    figname='SpecShapes'
    SetPaperFigPath
    print( fullfile(figdir ,figname ) , '-dpng' )
end

%%