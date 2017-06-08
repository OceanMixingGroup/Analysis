%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotResampT2.m
%
%
% 2 Mar 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot epsilon profile from true and one test

clear ; close all

testnum=9
minOT=10

useOrig=0
useGridded=1 ; dzgrid=10 ; dtgrid=2;

DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]

if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
    % load 'true' data
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end

load(fname)

em=nanmean(REsamp.eps,2);
figure(1);clf
agutwocolumn(0.8);wysiwyg


ht=semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)
hold on
semilogx(squeeze(em),REsamp.z,'.');axis ij
hold on
grid on
hs=semilogx(nanmean(squeeze(em),2),REsamp.z,'linewidth',2);axis ij
xlabel('<\epsilon>')
ylabel('Depth')
xlim([1e-10 10^(-6.5)])
%legend([ht hs],'true','samp','location','best')

% [bb]=bootstrap_profile(xx2.eps,100);
% semilogx(bb(:,2),xx2.z,'k','linewidth',2)
% semilogx(bb(:,[1 3]),xx2.z,'k-','linewidth',1)
%
clear bb
[a b c]=size(REsamp.eps)
e2=nan*ones(a,b*c);
for whz=1:a
    X=squeeze(REsamp.eps(whz,:,:));
   e2(whz,:)=X(:)'; 
end
%
%[bb]=bootstrap_profile(e2,100);
% h=semilogx(bb(:,2),REsamp.z)
% set(h,'Color',0.5*[1 1 1],'LineWidth',2)
% h=semilogx(bb(:,[1 3]),REsamp.z);
% set(h,'Color',0.5*[1 1 1],'LineStyle','-','LineWidth',1)

%% Plot profiles from a couple tests on same figure

clear ; close all

minOT=10

useOrig=0
useGridded=1 ; dzgrid=10 ; dtgrid=2;

clear testnum DatDir fname em REsamp
testnum=3
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
%fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end
load(fname)
em=nanmean(REsamp.eps,2);

figure(1);clf
agutwocolumn(0.8);wysiwyg

if useOrig==1
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
elseif useGridded==1
    % load 'true' data
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
end

idreal=isin(xx2.yday,[nanmin(REsamp.tgrid(:)) nanmax(REsamp.tgrid(:))]);
%semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)
semilogx(nanmean(xx2.eps(:,idreal),2),xx2.z,'k','linewidth',2)

hold on

%semilogx(squeeze(em),REsamp.z,'.');axis ij
grid on
semilogx(nanmean(squeeze(em),2),REsamp.z,'linewidth',2);axis ij

clear testnum DatDir fname em REsamp
testnum=1
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
%fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end

load(fname)
em=nanmean(REsamp.eps,2);
semilogx(nanmean(squeeze(em),2),REsamp.z,'linewidth',2);axis ij


clear testnum DatDir fname em REsamp
testnum=2
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
%fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end

load(fname)
em=nanmean(REsamp.eps,2);
semilogx(nanmean(squeeze(em),2),REsamp.z,'linewidth',2);axis ij

clear testnum DatDir fname em REsamp
testnum=9
DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
%fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end

load(fname)
em=nanmean(REsamp.eps,2);
semilogx(nanmean(squeeze(em),2),REsamp.z,'linewidth',2);axis ij

xlim([1e-10 10^(-7)])
xlabel('<\epsilon>','fontsize',16)
ylabel('Depth','fontsize',16)

legend('true','15cm/s','25cm/s','50cm/s')


fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['T2_EpsProfResamp'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-png')

%%


clear ; close all

testnum=3
minOT=10

useOrig=0
useGridded=1 ; dzgrid=10 ; dtgrid=2;

DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]

if useOrig==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
    % load 'true' data
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
elseif useGridded==1
    fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end

load(fname)
whc=55

cl=[-9 -5]

figure(1);clf
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,2);

tm=nanmean(xx2.T,2);

axes(ax(1))
ezpc(xx2.yday,xx2.z,log10(xx2.eps))
colorbar
caxis(cl)
cmap=jet;
colormap([ 0.8*[1 1 1] ; cmap])
hold on
contour(xx2.yday,xx2.z,xx2.T,tm(1:5:end),'k')
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')
xlim([18.5 21.2])

axes(ax(2))
ezpc(REsamp.tgrid(whc,:),REsamp.z,log10(REsamp.eps(:,:,whc)))
colorbar
caxis(cl)
cmap=jet;
colormap([ 0.8*[1 1 1] ; cmap])
xlabel('Yearday')

linkaxes(ax)

%%