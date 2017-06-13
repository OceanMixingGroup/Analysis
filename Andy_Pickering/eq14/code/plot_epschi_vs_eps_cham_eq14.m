%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot eps_chi vs eps for Chameleon data.
%
%
%~~~~~~~~~~~~~~~~~
% 6/13/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

eps_chi = cham.CHI .* cham.N2 / 2 / 0.2 ./ (cham.DTDZ.^2);

figure(1);clf
agutwocolumn(0.6)
wysiwyg

histogram2( log10(cham.EPSILON), log10(eps_chi), 'DisplayStyle','tile')
xlim([-10 -5])
ylim([-12 -5])
xvec = linspace(-10,-5,100);
hold on
plot(xvec,xvec,'k--')
xlabel('\epsilon','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)

eq14_patches_paths
figname = [project_short '_epschi_vs_eps']
print(fullfile(fig_dir, figname), '-dpng')

%%