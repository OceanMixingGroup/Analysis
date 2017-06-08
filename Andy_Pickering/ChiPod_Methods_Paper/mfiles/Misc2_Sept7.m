%%
%
% Misc2_Sept7.m
%
%
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/TestChiMethod/

% load EQ14 data
dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';

clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )

gam1=est_gam(cham);
%%
clear idb
idb=find(log10(cham.EPSILON)<-8.25);
cham.CHI(idb)=nan;
cham.EPSILON(idb)=nan;
gam1=est_gam(cham);
%clear cham

%

idb=find(gam1>1);
gam1(idb)=nan;
clear xbins ybins
%xbins=-12:0.1:-4;
xbins=-7:0.05:-2;
ybins=-4:0.05:1 ;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/

% Plot gamma vs N2
figure(1);clf
agutwocolumn(0.6)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(cham.N2(:)),0,log10(gam1(:)),0,0);
h=pcolor(xbins,ybins,hist);
set(h,'edgecolor','none')
hold on
 plot(xbins,xbins+2,'r--')
grid on
cmap=flipud(hot);
colormap([0.85*[1 1 1] ; cmap])
xlabel('log_{10} N^2','fontsize',16)
ylabel('log_{10} \Gamma','fontsize',16)
%%

figure(2);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(411)
ezpc(cham.castnumber,cham.depth,real(log10(cham.N2)))
colorbar
ylim([0 200])
SubplotLetterMW('log_{10}[N^2]')


ax2=subplot(412)
ezpc(cham.castnumber,cham.depth,real(log10(cham.DTDZ)))
colorbar
caxis([-3 1])
SubplotLetterMW('log_{10}[dT/dz]')

ax3=subplot(413)
ezpc(cham.castnumber,cham.depth,real(log10(cham.EPSILON)))
colorbar
SubplotLetterMW('log_{10}[\epsilon]')

ax4=subplot(414)
ezpc(cham.castnumber,cham.depth,gam1)
colorbar
caxis([0 0.1])
SubplotLetterMW('log_{10}[\Gamma]')

linkaxes([ax1 ax2 ax3 ax4])

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures'
figname='pcolor_n2_dtdz_eps_gam_eq14_chameleon'
print(fullfile(figdir,figname),'-dpng')

%%
ylim([0 200])
xlim([1600 2000])
shg
figname='pcolor_n2_dtdz_eps_gam_eq14_chameleon_zoom'
print(fullfile(figdir,figname),'-dpng')
%%
gam1(find(gam1>1))=nan;
figure(8);clf
agutwocolumn(0.75)
wysiwyg

plot(gam1(:),cham.P(:),'.')
hold on
plot(nanmean(gam1,2),nanmean(cham.P,2),'ko')
%xlim([0 2])
axis ij
xlabel('\Gamma')
ylabel('P')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures'
figname='Hist_gamma_vs_P_eq14_chameleon'
print(fullfile(figdir,figname),'-dpng')

%%
