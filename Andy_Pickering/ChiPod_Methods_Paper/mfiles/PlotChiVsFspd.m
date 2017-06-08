%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotChiVsFspd.m
%
% Plot chi from EQ14 CTD-chipod for different ranges of u to show there
% isn't a bias.
%
%---------------
% 05/26/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'

saveplots=1

uall=[];
chiall=[];
hb=waitbar(0)
for castnum=1:3100
    waitbar(castnum/3100,hb)
    fmax=7;
    try
        clear avg
        load( fullfile(datdir,['zsm10m_fmax' num2str(fmax) 'Hz_respcorr0_fc_99hz_gamma20'],['EQ14_' sprintf('%04d',castnum) 'avg.mat']) )
        %avg1=avg;clear avg5
        
        uall=[uall ; avg.fspd(:)];
        chiall=[chiall ; avg.chi1(:)];
        
    end
    
end
delete(hb)

%%

uall=abs(uall);


%%
figure(1);clf
loglog(uall,chiall,'.')
xlabel('U')
ylabel('\chi')
xlim([1e-1 10^(0.2)])
%%
figure(1);clf
loglog(avg.fspd,avg.chi1,'.')
grid on

%% try splitting up by U ranges

ig0=find( uall<0.3 );
ig1=find(uall>0.3 & uall<0.5 );
ig2=find(uall>0.5 & uall<0.75);
ig3=find(uall>0.75 & uall<1);
ig4=find(uall>1 );

Ds='stair'

figure(1);clf
h1=histogram(log10(chiall(ig1)),'Normalization','pdf','DisplayStyle',Ds)
hold on
h2=histogram(log10(chiall(ig2)),h1.BinEdges,'Normalization','pdf','DisplayStyle',Ds)
h3=histogram(log10(chiall(ig3)),h1.BinEdges,'Normalization','pdf','DisplayStyle',Ds)
%h3=histogram(log10(chiall(ig4)),h1.BinEdges,'Normalization','pdf','DisplayStyle',Ds)
h0=histogram(log10(chiall(ig0)),h1.BinEdges,'Normalization','pdf','DisplayStyle',Ds)
grid on
xlabel('log_{10} \chi')
legend([h0 h1 h2 h3],'u<0.3','0.3>u<0.5','0.5>u>0.75','u>0.75')
title('EQ14 All CTD chi-pod profiles')

%%

SetNotesFigDir
print(fullfile(NotesFigDir,'EQ14_chi_hist_diffU'),'-dpng')

%%

figure(1);clf
loglog(abs(uall(ig0)),chiall(ig0),'.')
hold on
loglog(abs(uall(ig1)),chiall(ig1),'.')
loglog(abs(uall(ig2)),chiall(ig2),'.')
loglog(abs(uall(ig3)),chiall(ig3),'.')
loglog(abs(uall(ig4)),chiall(ig4),'.')
grid on


%%

figure(2);clf
semilogy(abs(uall(ig0)),chiall(ig0),'.')
hold on
semilogy(abs(uall(ig1)),chiall(ig1),'.')
semilogy(abs(uall(ig2)),chiall(ig2),'.')
semilogy(abs(uall(ig3)),chiall(ig3),'.')
semilogy(abs(uall(ig4)),chiall(ig4),'.')
grid on
xlim([0 1.5])
xlabel('abs(U)','fontsize',16)
ylabel('log_{10}\chi','fontsize',16)

%%
SetNotesFigDir
print(fullfile(NotesFigDir,'eq14cham_chi_vs_U'),'-dpng')

%%


clear xbins ybins
%xbins=-12:0.1:-4;
xbins=-0.5:0.025:0.5;
ybins=-12:0.05:-4;

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
[hist,mn,mdn,md]=hist2d(xbins,ybins,log10(uall),0,log10(chiall),0,0);
%
figure(9);clf
agutwocolumn(0.5)
wysiwyg
h=pcolor(xbins,ybins,hist);
set(h,'edgecolor','none')
hold on
plot(xbins,xbins,'k--')
plot(xbins,xbins-1,'r--')
plot(xbins,xbins+1,'r--')
colorbar
grid on
cmap=flipud(hot);
colormap([0.85*[1 1 1] ; cmap])

%%