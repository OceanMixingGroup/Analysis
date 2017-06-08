%~~~~~~~~~~~~~~~~~~~~~
%
% PlotPctDataDepthLoop.m
%
% Make a figure showing the % of data that is thrown out as depth loops for
% each cast in EQ14 CTD-chipod, and the % of depth ranges that we still
% have good data for. 
%
%---------------
% 05/27/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%Flist=dir(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma5/','*T1*'))
Flist=dir(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20/','*downcast*T1*'))
Flist=dir(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20/','*upcast*T1*'))

pct_all=[];
dgd=[]
for icast=1:length(Flist)

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma5/',Flist(icast).name))
pct_all=[pct_all avg.pct_remove];
ig=find(avg.good_depth_bins>5);
dgd=[dgd length(ig)/length(avg.good_depth_bins)*100];
end

figure(1);clf
agutwocolumn(0.6)
wysiwyg

subplot(121)
histogram(pct_all,'BinEdges',0:10:100)
xlabel('% data thrown out per cast')
%xlim([0 50])

subplot(122)
histogram(dgd,'BinEdges',0:10:100)
xlabel('% 10m bins w/ >5 points')
%xlim([40 100])

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper'
print(fullfile(figdir,['EQ14_pctDloop']),'-dpng')

%%