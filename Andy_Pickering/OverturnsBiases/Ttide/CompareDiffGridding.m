%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareDiffGridding.m
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% compare overturns computed from raw and gridded data

clear ; close all

minOT=10

% ~~ compare same tgrid, different zgrids
clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X1=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=15  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X2=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=20  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X3=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=30  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X4=xx2;clear xx2
%~~
%

figure(1);clf
h1=semilogx(nanmean(X1.eps,2),X1.z,'linewidth',2)
hold on
h2=semilogx(nanmean(X2.eps,2),X2.z,'linewidth',2)
h3=semilogx(nanmean(X3.eps,2),X3.z,'linewidth',2)
h4=semilogx(nanmean(X4.eps,2),X4.z,'linewidth',2)

load /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_10
h0=semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)

legend

axis ij
grid on
xlabel('<\epsilon>','fontsize',16)
xlim([1e-10 10^(-6.8)])
ylabel('Depth','fontsize',16)
legend([h0 h1 h2 h3 h4],'Orig','10m','15m','20m','30m')


fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['T2_EpsProfile_diffdz'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-png')

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

ax = MySubplot(0.1, 0.03, 0.07, 0.06, 0.1, 0.003, 1,5);

t_avg=1/24
[e0,t0]=SimpleBoxCar(nansum(xx2.eps),t_avg,xx2.yday);
[e1,t1]=SimpleBoxCar(nansum(X1.eps),t_avg,X1.yday);
[e2,t2]=SimpleBoxCar(nansum(X2.eps),t_avg,X2.yday);
[e3,t3]=SimpleBoxCar(nansum(X3.eps),t_avg,X3.yday);
[e4,t4]=SimpleBoxCar(nansum(X4.eps),t_avg,X4.yday);

axes(ax(1))
%semilogy(xx2.yday,nansum(xx2.eps),'k')
semilogy(t0,e0,'k')
hold on
semilogy(t1,e1)
semilogy(t2,e2)
semilogy(t3,e3)
semilogy(t4,e4)
% semilogy(X1.yday,nansum(X1.eps))
% semilogy(X2.yday,nansum(X2.eps))
% semilogy(X3.yday,nansum(X3.eps))
% semilogy(X4.yday,nansum(X4.eps))
%ylim([1e-10 1e-5])
axis tight
grid on
cb=colorbar;killcolorbar(cb)

axes(ax(2))
ezpc(X1.yday,X1.z,log10(X1.eps))
colorbar
caxis([-9 -5])

axes(ax(3))
ezpc(X2.yday,X2.z,log10(X2.eps))
colorbar
caxis([-9 -5])

axes(ax(4))
ezpc(X3.yday,X3.z,log10(X3.eps))
colorbar
caxis([-9 -5])

axes(ax(5))
ezpc(X4.yday,X4.z,log10(X4.eps))
colorbar
caxis([-9 -5])
cmap=jet;
colormap([ 0.8*[1 1 1] ; cmap])

linkaxes(ax,'x')

%%



h=figure(1);clf
agutwocolumn(0.8)
wysiwyg

ax = MySubplot(0.1, 0.03, 0.07, 0.06, 0.1, 0.07, 2,2);

axes(ax(1))
histogram(X1.Lot_each(:),20)
hold on
histogram(X2.Lot_each(:),20)
histogram(X3.Lot_each(:),20)
histogram(X4.Lot_each(:),20)
xlabel('L_{ot}')

axes(ax(2))
histogram(X1.Lt_each(:),20)
hold on
histogram(X2.Lt_each(:),20)
histogram(X3.Lt_each(:),20)
histogram(X4.Lt_each(:),20)
xlabel('L_{ot}')

%%
%% ~~ compare same zgrid, different tgrids
clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=1
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X1=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X2=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=5
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X3=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=10
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
X4=xx2;clear xx2
%~~

%%
%%

figure(2);clf

Nm='pdf'
Ds='stair'

whvar='Otnsq_each'
whvar='Lot_each'
whvar='Lt_each'

figure(2);clf
h1=histogram(X1.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
hold on
h2=histogram(X2.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
h3=histogram(X3.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
h4=histogram(X4.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
h0=histogram(xx2.(whvar)(:),40,'Normalization',Nm,'DisplayStyle',Ds,'EdgeColor','k')
grid on
% freqline(nanmean(X1.(whvar)(:)))
% freqline(nanmean(X2.(whvar)(:)))
% freqline(nanmean(X3.(whvar)(:)))
xlabel(whvar,'interpreter','none')

legend([h0 h1 h2 h3 h4],'Orig','1','2','3','4','location','best')
%%
whvar='eps_each'
Ds='bar'
Ds='stair'
figure(3);clf
h1=histogram(log10(X1.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
hold on
h2=histogram(log10(X2.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
h3=histogram(log10(X3.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
h4=histogram(log10(X4.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
h0=histogram(log10(xx2.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds,'EdgeColor','k')
grid on
%
% freqline(log10(nanmean(X1.(whvar)(:))))
% freqline(log10(nanmean(X2.(whvar)(:))))
% freqline(log10(nanmean(X3.(whvar)(:))))
xlabel('log_{10} \epsilon')
legend([h0 h1 h2 h3 h4],'Orig','1','2','3','4','location','best')
%%

clear ; close all

figure(1);clf

minOT=10

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)

hold on
%
% clear dzgrid dtgrid DatDir fname xx2
% dzgrid=10  ;  dtgrid=5
% DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
% fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
% load(fullfile(DatDir,fname))
% semilogx(nanmean(xx2.eps,2),xx2.z)

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=5
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)

load /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_10
semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)

%load /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_10_gridded
%semilogx(nanmean(xx2.eps,2),xx2.z,'k')

%legend('raw','grid')
legend

axis ij
grid on
xlabel('<\epsilon>')
ylabel('Depth')

%%

clear dzgrid dtgrid DatDir fname xx2
dzgrid=10  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)
X1=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=15  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)
X2=xx2;clear xx2

clear dzgrid dtgrid DatDir fname xx2
dzgrid=20  ;  dtgrid=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2']
fname=['T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']
load(fullfile(DatDir,fname))
semilogx(nanmean(xx2.eps,2),xx2.z)
X3=xx2;clear xx2

Nm='pdf'
Ds='stair'

whvar='Otnsq_each'
whvar='Lot_each'
%whvar='Lt_each'
figure(2);clf
histogram(X1.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(X2.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
histogram(X3.(whvar)(:),'Normalization',Nm,'DisplayStyle',Ds)
grid on
freqline(nanmean(X1.(whvar)(:)))
freqline(nanmean(X2.(whvar)(:)))
freqline(nanmean(X3.(whvar)(:)))
xlabel(whvar,'interpreter','none')
%%

whvar='eps_each'

figure(2);clf
histogram(log10(X1.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(log10(X2.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
histogram(log10(X3.(whvar)(:)),'Normalization',Nm,'DisplayStyle',Ds)
grid on

freqline(log10(nanmean(X1.(whvar)(:))))
freqline(log10(nanmean(X2.(whvar)(:))))
freqline(log10(nanmean(X3.(whvar)(:))))
xlabel('log_{10} \epsilon')
%
%
%% Compare resampled OT for gridded and non gridded

clear ; close all


load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test1/T2_Test1_minOT_10_AllCases_Gridded_dz_10m_dt_2min')
em=squeeze(nanmean(REsamp.eps,2));

figure(1);clf
semilogx(nanmean(em,2),REsamp.z)
axis ij
hold on

clear REsamp em
load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test1/T2_Test1_minOT_10_AllCases')
em=squeeze(nanmean(REsamp.eps,2));
semilogx(nanmean(em,2),REsamp.z)
%%