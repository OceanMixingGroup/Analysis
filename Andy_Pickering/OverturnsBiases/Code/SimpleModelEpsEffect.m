%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% SimpleModelEpsEffect.m
%
% Simple model for effect on epsilon of symmetrically under/over estimating
% overturn sizes. 
%
% Original 22 Nov
% 17 Dec - Rename, comment
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Np=50
Lt=linspace(2,100,Np);
alpha=1.4;

eps_real=Lt.^2;
eps_samp=  (  (alpha.*Lt).^2 + (Lt./alpha).^2  ) ./2;

figure(1);clf
%plot(eps_real,eps_real,eps_real,eps_samp)
plot(Lt,eps_real,Lt,eps_samp)
%loglog(Lt,eps_real,Lt,eps_samp)
legend('real','samp','location','best')
xlabel('L_T')
ylabel('\epsilon')

figure(2);clf
%plot(Lt,(eps_samp-eps_real)./eps_real)
plot(Lt,(eps_samp-eps_real))

%%