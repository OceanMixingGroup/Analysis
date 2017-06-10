%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_chi_diff_respcorr_EQ08.m
%
% This script looks at effect of doing a response correction parameter in chi-pod calculations.
% I ran the chipod method on chameleon FP07 data from eq08 for different
% values , and here plot some comparisons.
%
%
%-----------------
% 03/17/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% main folder where data is
datadir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/'

% folder to save figures to
figdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures'

saveplot=1

% use same parameters for all cases except zsmooth
Params.gamma=0.2
Params.fmax=20
Params.z_smooth=10;

% load data using different 'respcorr' values

clear fmax fname C
Params.resp_corr=0
fname=['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100) '.mat']
load(fullfile(datadir,fname))
C1=C;clear C

clear fmax fname C
Params.resp_corr=1
fname=['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100) '.mat']
load(fullfile(datadir,fname))
C3=C;clear C

%% 1st Plot some histograms of chi and epsilon for different zsmooth values to see how they vary

Ds='stair';

figure(1);clf

subplot(211)
h1=histogram(log10(C1.chi(:)),'DisplayStyle',Ds);
hold on
h3=histogram(log10(C3.chi(:)),h1.BinEdges,'DisplayStyle',Ds);
xlim([-13 -2])
grid on
xlabel('log_{10}\chi')
legend([ h1  h3 ],'respcorr=0','respcorr=1','location','best')
title('EQ08 \chi pod method applied to chameleon')

subplot(212)
h1=histogram(log10(C1.eps(:)),'DisplayStyle',Ds);
hold on
h3=histogram(log10(C3.eps(:)),h1.BinEdges,'DisplayStyle',Ds);
xlim([-13 -2])
grid on
xlabel('log_{10}\epsilon')
legend([ h1  h3 ],'respcorr=0','respcorr=1','location','best')

if saveplot==1
    figname=['hist_chi_eps_diff_respcorr']
    print(fullfile(figdir,figname),'-dpng')
end


%% make scatter plot of variables from 2 of the cases

whvar='dTdz'
whvar='chi'

figure(1);clf
loglog(C1.(whvar)(:),C3.(whvar)(:),'.')
grid on
xvec=linspace(1e-12,1e0,100);
hold on
loglog(xvec,xvec,'k--')

%% average chameleon in same bins so we can make scatter plots etc.

% load Chameleon data
load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/sum/eq14_sum_clean_new_cstar.mat')

% parameters for binning
clear zmin dz zmax tbin zbin sbin
zmin=0;
zmax=300;
dz=5;
minobs=1;
zbins=zmin:dz:zmax ;

clear chi2 cham_eps N22 dTdz2
cham_chi =nan*ones(length(zbins),length(cham.castnumber));
cham_eps =cham_chi;
cham_N2  =cham_chi;
cham_dTdz=cham_chi;

for icast=1:length(C1.castnumber)
    ig=find(cham.castnumber==C1.castnumber(icast));
    %     [chi2(:,icast)  zbin Nobs] = binprofile( cham.CHI(:,icast)     , cham.P(:,icast), zmin, dz, zmax,minobs) ;
    %     [eps2(:,icast)  zbin Nobs] = binprofile( cham.EPSILON(:,icast) , cham.P(:,icast), zmin, dz, zmax,minobs) ;
    %     [N22(:,icast)   zbin Nobs] = binprofile( cham.N2(:,icast)      , cham.P(:,icast), zmin, dz, zmax,minobs) ;
    %     [dTdz2(:,icast) zbin Nobs] = binprofile( cham.DTDZ(:,icast)    , cham.P(:,icast), zmin, dz, zmax,minobs) ;
    
    [cham_chi(:,icast)  zbin Nobs] = binprofile( cham.CHI(:,ig)     , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_eps(:,icast)  zbin Nobs] = binprofile( cham.EPSILON(:,ig) , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_N2(:,icast)   zbin Nobs] = binprofile( cham.N2(:,ig)      , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_dTdz(:,icast) zbin Nobs] = binprofile( cham.DTDZ(:,ig)    , cham.P(:,ig), zmin, dz, zmax,minobs) ;
end



%% 2D hist of chameleon vs chipod method - chi



xbins=-12:0.05:-4;
ybins=xbins;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/

figure(1);clf
agutwocolumn(1)
wysiwyg

xl=[-12 -4]
cl=[0 130]

cmap=flipud(hot);
cmap=[ 0.9*[1 1 1] ; cmap ]

%
%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_chi(:)),0,log10(C1.chi(:)),0,0);

subplot(211)
%subplot(222)
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
title(['No response correction'])
xlabel('cham : log_{10}\chi ','fontsize',16)
 ylabel('\chi pod method : log_{10}\chi ','fontsize',16)
caxis(cl)


%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_chi(:)),0,log10(C3.chi(:)),0,0);

subplot(212)
%subplot(224)
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
title([' response correction'])
xlabel('cham : log_{10}\chi ','fontsize',16)
 ylabel('\chi pod method : log_{10}\chi ','fontsize',16)
caxis(cl)

if saveplot==1
    figname=['2Dhist_chi_diff_respcorr']
    print(fullfile(figdir,figname),'-dpng')
end

%
%% 2D hist of chameleon vs chipod method - eps

xbins=-12:0.05:-4;
ybins=xbins;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/

%
figure(1);clf
agutwocolumn(1)
wysiwyg

xl=[-12 -4]
cl=[10 80]

cmap=flipud(hot);
cmap=[ 0.9*[1 1 1] ; cmap ]


%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_eps(:)),0,log10(C1.eps(:)),0,0);

subplot(211)
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
xlabel('cham : log_{10}\epsilon ','fontsize',16)
title(['No response correction'])
caxis(cl)
 ylabel('\chi pod method : log_{10}\epsilon ','fontsize',16)
 

%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_eps(:)),0,log10(C3.eps(:)),0,0);

subplot(212)
%subplot(224)
h=pcolor(xbins,ybins,hist)
%hold on
%contour(xbins,ybins,hist,[20 100],'k')
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
xlabel('cham : log_{10}\epsilon ','fontsize',16)
title([' response correction'])
caxis(cl)
 ylabel('\chi pod method : log_{10}\epsilon ','fontsize',16)
 
if saveplot==1
    figname=['2Dhist_eps_diff_respcorr']
    print(fullfile(figdir,figname),'-dpng')
end


%%