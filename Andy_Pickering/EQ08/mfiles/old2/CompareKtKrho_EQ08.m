%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareKtKrho_EQ08.m
%
% Compute K_T and K_rho from EQ08 chameleon data to test the chipod assumption
% that they are equal.
%
% Modified from CompareKtKrho.m
%
%----------------
% 03/17/16 - AP - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load Chameleon data
load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')

% compute KT
KT=0.5 * cham.CHI ./ cham.DTDZ.^2 ;

% compute Krho
gamma=0.2
Krho=gamma * cham.EPSILON ./ cham.N2;

%% 
xl=[1e-12 1e2]
figure(1);clf
agutwocolumn(0.6)
wysiwyg
loglog(Krho(:),KT(:),'.')
xlabel('Krho')
ylabel('KT')
xvec=linspace(1e-10,1e1,100);
hold on
xlim(xl)
ylim(xl)
loglog(xvec,xvec,'k--')
loglog(xvec,0.8*xvec,'r--')
grid on

%%
xbins=-10:0.025:1;
ybins=xbins;

%figure(2);clf
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(Krho(:)),0,log10(KT(:)),0,0);
%[hist,mn,mdn,md]=hist2d(xbins,ybins,chi.N2(:),0,N22(:),0,0);
%%

figure(2);clf
agutwocolumn(0.6)
wysiwyg

h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
colorbar
grid on
cmap=flipud(hot);   
colormap([ 0.75*[1 1 1] ;cmap ])
caxis([1 90])
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlabel('log_{10}K_{\rho}','fontsize',16)
ylabel('log_{10}K_T','fontsize',16)
title('EQ08 Chameleon Data','fontsize',16)

%%

print('/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures/Cham_KtVsKrho_2Dhist','-dpng')

%%

figure(1);clf
histogram(log10(KT./Krho),50)
xlabel('log_{10} [K_T / K_{\rho}]','fontsize',16)
grid on

%%

print('/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures/Cham_Kt_Krho_hist','-dpng')

%%