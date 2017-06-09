%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_ChiProc_To_Cham_eq08_AP.m
%
% Make plots of chipod-method results vs actual for EQ08 chameleon data.
%
% Modified from Plot_chiProc_vs_actual_EQ14.m
%
%---------------------
% 03/16/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=1

Params.fmax=10;
Params.z_smooth=10;
Params.gamma=0.2;
Params.resp_corr=1;
Params.nfft=128

%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chi_all_' num2str(z_smooth) 'm.mat'])
%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chi_all_' num2str(Params.z_smooth) 'm_fmax_' num2str(Params.fmax) '.mat'])

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/',...
   ['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100) '_nfft' num2str(Params.nfft)]));


load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')

figdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Figures'

%% Pcolor N2 and dTdz comparison

xl=[0 3100]
yl=[0 250]

cl=[-5 -3]

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(411);
ezpc(C.castnumber,C.p,real(log10(C.N2)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('N2 chipod method')
ylabel('pres. [db]','fontsize',16)

ax2=subplot(412);
ezpc(cham.castnumber,cham.P,real(log10(cham.N2)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('N2 chameleon')
ylabel('pres. [db]','fontsize',16)

cl=[-3 0]
ax3=subplot(413);
ezpc(C.castnumber,C.p,real(log10(C.dTdz)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('dTdz chipod method')
ylabel('pres. [db]','fontsize',16)

ax4=subplot(414);
ezpc(cham.castnumber,cham.P,real(log10(cham.DTDZ)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('dTdz chameleon')
ylabel('pres. [db]','fontsize',16)
xlabel('castnumber','fontsize',16)

linkaxes([ax1 ax2 ax3 ax4])

if saveplot==1
    
    print( fullfile( figdir,['EQ14_chiVsTrue_N2_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)] ) , '-dpng' )
end


%% Pcolor chi comparison

%saveplot=0

xl=[0 2650];
yl=[0 200];
cl=[-10 -5];

figure(2);clf
agutwocolumn(0.8)
wysiwyg

ax1=subplot(211);
ezpc(C.castnumber,C.p,log10(C.chi))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('chipod method \chi','fontsize',16);
ylabel('Pres [db]','fontsize',16)

ax2=subplot(212);
ezpc(cham.castnumber,cham.P,log10(cham.CHI))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('chameleon \chi','fontsize',16);
xlabel('castnumber','fontsize',16)
ylabel('Pres [db]','fontsize',16)

linkaxes([ax1 ax2])

if saveplot==1
    figname=['EQ08_chiVsTrue_chi_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname) , '-dpng' )
     % save copy to paper directory
    figdir2='/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper'
    print( fullfile( figdir2,figname) , '-dpng' )
end



%% Pcolor eps comparison

xl=[0 3100];
yl=[0 250];
cl=[-11 -4];

figure(3);clf
agutwocolumn(0.8)
wysiwyg

ax1=subplot(211);
ezpc(C.castnumber,C.p,real(log10(C.eps)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('chipod method \epsilon','fontsize',16);
ylabel('Pres [db]','fontsize',16)

ax2=subplot(212);
ezpc(cham.castnumber,cham.P,real(log10(cham.EPSILON1)))
colorbar
xlim(xl)
ylim(yl)
caxis(cl)
title('chameleon \epsilon','fontsize',16)
xlabel('castnumber','fontsize',16)
ylabel('Pres [db]','fontsize',16)

linkaxes([ax1 ax2])

if saveplot==1
    figname=['EQ08_chiVsTrue_eps_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end


%% Plot histograms of chi and epsilon from both methods

%saveplot=0

Nm='pdf';
figure(4);clf
agutwocolumn(0.7)
wysiwyg
m=1,n=2;

subplot(m,n,1)
h1=histogram(log10(C.chi(:)),'DisplayStyle','Stair','Normalization',Nm);
hold on
histogram(log10(cham.CHI(:)),'BinEdges',h1.BinEdges,'DisplayStyle','Stair','Normalization',Nm)
xlim([-13 0])
legend('\chi pod','cham','location','best')
grid on
xlabel('log_{10} \chi','fontsize',16)
ylabel(Nm)

subplot(m,n,2)
h1=histogram(log10(C.eps(:)),'DisplayStyle','Stair','Normalization',Nm);
hold on
histogram(log10(cham.EPSILON(:)),'BinEdges',h1.BinEdges,'DisplayStyle','Stair','Normalization',Nm)
xlim([-13 -2])
legend('\chi pod','cham','location','best')
grid on
xlabel('log_{10} \epsilon','fontsize',16)

if saveplot==1
   figname=['EQ08_chiVsTrue_Hist_chieps_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end


%% histogram of N2 and dTdz from both

Nm='pdf';

figure(5);clf
agutwocolumn(0.5)
wysiwyg

ax1=subplot(121)
histogram(log10(C.N2(:)),'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(real(log10(cham.N2(:))),'DisplayStyle','Stair','Normalization',Nm)
xlim([-8 0])
legend('\chi pod','cham','location','best')
xlabel('log_{10}N^2')
ylabel(Nm)
grid on

ax2=subplot(122)
histogram(log10(C.dTdz(:)),'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(real(log10(cham.DTDZ(:))),'DisplayStyle','Stair','Normalization',Nm)
xlim([-6 1])
legend('\chi pod','cham','location','best')
xlabel('log_{10}dT/dz')
ylabel(Nm)
grid on

%%  Plot chi/eps vs dTdz etc for chipod method - 2D hist

xbins=-12:0.05:1;
ybins=xbins;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(C.dTdz(:)),0,log10(C.chi(:)),0,0);

xl=[-4 0]

figure(6);clf
agutwocolumn(1)
wysiwyg
m=3,n=1;

ax1=subplot(m,n,1);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')
xlim(xl)
ylim([-12 0])
shading flat
colorbar
cmap=flipud(hot);
cmap=[0.75*[1 1 1] ; cmap]
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\chi','fontsize',16)
title('\chi pod method','fontsize',16)
%
%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(C.dTdz(:)),0,log10(C.eps(:)),0,0);

ax2=subplot(m,n,2);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')
xlim(xl)
ylim([-12 0])
shading flat
colorbar
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)
title(['zsmooth=' num2str(Params.z_smooth) 'm, fmax=' num2str(Params.fmax) 'Hz , \Gamma =' num2str(Params.gamma)])

%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(C.dTdz(:)),0,log10(C.KT(:)),0,0);

ax3=subplot(m,n,3);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')

xlim(xl)
ylim([-12 0])
shading flat
colorbar
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}K_T','fontsize',16)

linkaxes([ax1 ax2 ax3],'x')

if saveplot==1
   figname=['EQ08_chimeth_vsdTdz_2Dhist_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end


%%  Plot chi/eps vs dTdz etc for chameleon - 2D hist

xbins=-12:0.05:1;
ybins=xbins;
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham.DTDZ(:)),0,log10(cham.CHI(:)),0,0);

xl=[-4 0]

figure(7);clf
agutwocolumn(1)
wysiwyg
%xl=[1e-6 1e0];
m=3,n=1;

ax1=subplot(m,n,1);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')
xlim(xl)
ylim([-12 -2])
shading flat
colorbar
cmap=flipud(hot);
cmap=[0.85*[1 1 1] ; cmap]
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\chi','fontsize',16)
title('chameleon','fontsize',16)
%
%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham.DTDZ(:)),0,log10(cham.EPSILON(:)),0,0);

ax2=subplot(m,n,2);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')
xlim(xl)
ylim([-12 0])
shading flat
colorbar
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)

%~~
clear hist mn mdn md
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham.DTDZ(:)),0,log10(0.5*cham.CHI(:)./(cham.DTDZ(:).^2)),0,0);

ax3=subplot(m,n,3);
h=pcolor(xbins,ybins,hist)
hold on
plot(-4:0,(-4:0)-4,'k--')
plot(-4:0,-(-4:0)-8,'k--')

xlim(xl)
ylim([-12 0])
shading flat
colorbar
colormap(cmap)
grid on
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}K_T','fontsize',16)

linkaxes([ax1 ax2 ax3],'x')

if saveplot==1
    figname=['EQ08_cham_vsdTdz_2Dhist']
    print( fullfile( figdir,figname ) , '-dpng' )
end


%% Make a scatter plot comparing the two methods; 1st average chameleon data in same bins

% * 03/16/16 - AP - cham castnumbers not exactly equal those of chi proc...

% parameters for binning
clear zmin dz zmax tbin zbin sbin
zmin=0;
dz=5;
zmax=300;
minobs=1;
zbins=zmin:dz:zmax ;

clear chi2 eps2 N22 dTdz2
%chi2=nan*ones(length(zbins),length(cham.castnumber));
chi2=nan*ones(length(zbins),length(C.castnumber));
eps2=chi2;
N22=chi2;
dTdz2=chi2;
%KT2=chi2;

for icast=1:length(C.castnumber)
    ig=find(cham.castnumber==C.castnumber(icast));
%     [chi2(:,icast)  zbin Nobs] = binprofile( cham.CHI(:,icast)     , cham.P(:,icast), zmin, dz, zmax,minobs) ;
%     [eps2(:,icast)  zbin Nobs] = binprofile( cham.EPSILON(:,icast) , cham.P(:,icast), zmin, dz, zmax,minobs) ;
%     [N22(:,icast)   zbin Nobs] = binprofile( cham.N2(:,icast)      , cham.P(:,icast), zmin, dz, zmax,minobs) ;
%     [dTdz2(:,icast) zbin Nobs] = binprofile( cham.DTDZ(:,icast)    , cham.P(:,icast), zmin, dz, zmax,minobs) ;

    [chi2(:,icast)  zbin Nobs] = binprofile( cham.CHI(:,ig)     , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [eps2(:,icast)  zbin Nobs] = binprofile( cham.EPSILON(:,ig) , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [N22(:,icast)   zbin Nobs] = binprofile( cham.N2(:,ig)      , cham.P(:,ig), zmin, dz, zmax,minobs) ;
    [dTdz2(:,icast) zbin Nobs] = binprofile( cham.DTDZ(:,ig)    , cham.P(:,ig), zmin, dz, zmax,minobs) ;
end


%% 2D hist of N2 from both methods

xbins=-6:0.05:-2.5;
ybins=xbins;
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(C.N2(:)),0,log10(N22(:)),0,0);
%[hist,mn,mdn,md]=hist2d(xbins,ybins,C.N2(:),0,N22(:),0,0);


figure(8);clf
agutwocolumn(0.6)
wysiwyg
%grid on
h=pcolor(xbins,ybins,hist)
%h=contourf(xbins,ybins,hist,30,'edgecolor','none')
set(h,'edgecolor','none')
colorbar
grid on
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
caxis([10 700])
hold on
plot(xbins,xbins,'k--')
xlabel('\chi pod method : log_{10}\N^2 ','fontsize',16)
ylabel('chameleon : log_{10}\N^2 ','fontsize',16)
title('2D hist of N^2 from both methods','fontsize',16)

if saveplot==1
    figname=['EQ08_2dhist_N2_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end

%% Plot histogram of chipod chi / cham chj

figure(1);clf

histogram(log10(C.chi(:)./chi2(:)))
xlim(3*[-1 1])
xlabel(['log_{10} \chi_{\chi}/\chi_{\epsilon}'])
grid on

%% 2D hist of chi from both methods

%xbins=-12:0.1:-4;
xbins=-12:0.05:-4;

ybins=xbins;
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(chi2(:)),0,log10(C.chi(:)),0,0);
%
figure(9);clf
agutwocolumn(0.6)
wysiwyg
h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
colorbar
grid on
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
caxis([0 110])

hold on
plot(xbins,xbins,'k--')
ylabel('\chi pod method : log_{10}\chi ','fontsize',16)
xlabel('chameleon : log_{10}\chi ','fontsize',16)
title('2D hist of \chi from both methods','fontsize',16)

if saveplot==1
    figname=['EQ08_2dhist_chi_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
     % also print a copy to paper folder
    print( fullfile('/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper' ,figname ) , '-dpng' )
end

%% 2D hist of epsilon from both methods

xbins=-12:0.05:-4;
ybins=xbins;
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(eps2(:)),0,log10(C.eps(:)),0,0);
%
figure(10);clf
agutwocolumn(0.6)
wysiwyg

h=pcolor(xbins,ybins,hist)
set(h,'edgecolor','none')
colorbar
grid on
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
caxis([0 90])
hold on
plot(xbins,xbins,'k--')
ylabel('\chi pod method : log_{10}\epsilon ','fontsize',16)
xlabel('chameleon : log_{10}\epsilon ','fontsize',16)
title('2D hist of \epsilon from both methods','fontsize',16)


if saveplot==1
    figname=['EQ08_2dhist_eps_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end



%% Plot histograms of chi and eps from both methods

figure(11);clf
agutwocolumn(0.6)
wysiwyg

subplot(211)
h1=histogram(log10(C.chi(:)),'DisplayStyle','stair')
hold on
h2=histogram(log10(chi2(:)),'BinEdges',h1.BinEdges,'DisplayStyle','stair')
xlim([-13 0])
ylabel('count','fontsize',16)
xlabel('log_{10}\chi','fontsize',16)
grid on
legend([h1 h2],'chipod','cham','location','best')

subplot(212)
h1=histogram(log10(C.eps(:)),'DisplayStyle','stair')
hold on
h2=histogram(log10(eps2(:)),'BinEdges',h1.BinEdges,'DisplayStyle','stair')
xlim([-13 -2])
ylabel('count','fontsize',16)
xlabel('log_{10}\epsilon','fontsize',16)
grid on
legend([h1 h2],'chipod','cham','location','best')


if saveplot==1
    figname=['EQ08_Hist_gridded_chi_eps_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100)]
    print( fullfile( figdir,figname ) , '-dpng' )
end

%%
