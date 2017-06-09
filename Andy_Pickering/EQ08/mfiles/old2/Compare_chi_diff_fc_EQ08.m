%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_chi_diff_fc_EQ08.m
%
% This script looks at effect of 'fc' parameter in chi-pod calculations.
% I ran the chipod method on chameleon FP07 data from eq08 for different
% values of fc, and here plot some comparisons.
%
%
%-----------------
% 04/28/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% main folder where data is
datadir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/'

figdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures'

saveplot=0

% use same z_smooth for all cases
Params.z_smooth=10;
Params.gamma=0.2
Params.fmax=15;
Params.resp_corr=1
Params.nfft=128;

% load data using different 'fc' values
clear fmax fname C
Params.fc=10
fname=['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
load(fullfile(datadir,fname))
C1=C;clear C

%
clear fmax fname C
Params.fc=15
fname=['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
load(fullfile(datadir,fname))
C2=C;clear C

% load data with no correction to compare also
% clear fmax fname C
% Params.resp_corr=0
% Params.fc=99;
% fname=['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]
% load(fullfile(datadir,fname))
% C3=C;clear C


%% 1st Plot some histograms of chi and epsilon for different fmax values to see how they vary

Ds='stair';

figure(1);clf

subplot(211)

h1=histogram(log10(C1.chi(:)),'DisplayStyle',Ds);
hold on
h2=histogram(log10(C2.chi(:)),h1.BinEdges,'DisplayStyle',Ds);
xlim([-13 -2])
grid on
xlabel('log_{10}\chi','fontsize',16)
legend([h1 h2],'fc=10hz','fc=15hz','location','best')

subplot(212)
h1=histogram(log10(C1.eps(:)),'DisplayStyle',Ds);
hold on
h2=histogram(log10(C2.eps(:)),h1.BinEdges,'DisplayStyle',Ds);
xlim([-13 -2])
grid on
xlabel('log_{10}\epsilon','fontsize',16)
legend([h1 h2],'fc=10hz','fc=15hz','location','best')

%% make scatter plot of variables from 2 of the cases

whvar='dTdz'
whvar='chi'

figure(1);clf
loglog(C1.(whvar)(:),C2.(whvar)(:),'.')
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
    [cham_chi(:,icast)  zbin Nobs] = binprofile( cham.CHI(:,ig)     , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_eps(:,icast)  zbin Nobs] = binprofile( cham.EPSILON(:,ig) , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_N2(:,icast)   zbin Nobs] = binprofile( cham.N2(:,ig)      , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [cham_dTdz(:,icast) zbin Nobs] = binprofile( cham.DTDZ(:,ig)    , cham.P(:,ig), zmin, dz, zmax,minobs) ;
end


%% Make a Pcolor plot

cl=[-11 -3]

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(311);
ezpc(C1.castnumber,zbin,log10(cham_chi))
colorbar
caxis(cl)
ylim([0 200])
title('EQ Chameleon - log_{10}\chi')
%xlabel('castnumber','fontsize',16)

ax2=subplot(312);
ezpc(C1.castnumber,zbin,log10(C1.chi))
colorbar
caxis(cl)
ylim([0 200])
title(['fc=' num2str(C1.Params.fc) 'hz'])
%xlabel('castnumber','fontsize',16)
ylabel('pres [db]','fontsize',16)

ax3=subplot(313);
ezpc(C2.castnumber,zbin,log10(C2.chi))
colorbar
caxis(cl)
ylim([0 200])
title(['fc=' num2str(C2.Params.fc) 'hz'])
xlabel('castnumber','fontsize',16)

linkaxes([ax1 ax2 ax3])

colormap(parula)

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
cmap=[ 0.9*[1 1 1] ; cmap ];


%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_chi(:)),0,log10(C1.chi(:)),0,0);

subplot(211)
h=pcolor(xbins,ybins,hist);
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
title(['fc=' num2str(C1.Params.fc) 'hz'])
xlabel('cham : log_{10}\chi ','fontsize',16)
 ylabel('\chi pod method : log_{10}\chi ','fontsize',16)
caxis(cl)


%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham_chi(:)),0,log10(C2.chi(:)),0,0);

subplot(212)
h=pcolor(xbins,ybins,hist);
set(h,'edgecolor','none')
grid on
colormap(cmap)
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
xlim(xl);ylim(xl)
title(['fc=' num2str(C2.Params.fc) 'hz'])
xlabel('cham : log_{10}\chi ','fontsize',16)
ylabel('\chi pod method : log_{10}\chi ','fontsize',16)
caxis(cl)

if saveplot==1
    figname=['2Dhist_chi_' num2str(Params.z_smooth) 'mzsmooth_diff_fc']
    print(fullfile(figdir,figname),'-dpng')
    SetNotesFigDir
    print(fullfile(NotesFigDir,figname),'-dpng')
end

%


%%