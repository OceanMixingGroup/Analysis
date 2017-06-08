%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotHist_fspdErr.m
%
% * Makes plot for chipod methods paper *
%
% Formerly part of Plot_chi_diff_fspd.m
%
% **** run processing againf ro no resp corr
%
%---------------
% 04/18/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=1

perr1=[]
perr2=[]

for castnum=1:50
    
    try
        
        clear avg avg0
        datdir0='/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/';
        load(fullfile(datdir0,'zsm10m_fmax10Hz_respcorr0_fc_99hz_gamma20',['avg_cast' sprintf('%02d',castnum) '_upcast_SN1001_T1.mat']))
        avg0=avg;clear avg
        
        datdir='/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/avg/proc_varyfspd';
        
        % +.1 m/s
        whf=1;
        clear avg1 avg
        load(fullfile(datdir,['fvar_' num2str(whf)],['avg_cast' sprintf('%02d',castnum) '_upcast_SN1001_T1.mat']))
        avg1=avg;clear avg
        
        % +1 m/s
        whf=2;
        clear avg2 avg
        load(fullfile(datdir,['fvar_' num2str(whf)],['avg_cast' sprintf('%02d',castnum) '_upcast_SN1001_T1.mat']))
        avg2=avg;clear avg
        
        p1=100*(avg1.chi1-avg0.chi1)./avg0.chi1;
        p2=100*(avg2.chi1-avg0.chi1)./avg0.chi1;
        
        perr1=[perr1 ; p1(:) ];
        perr2=[perr2 ; p2(:) ];
        
    end % try
    
end % castnum
%%

perr1(find(abs(perr1)>200))=nan;
perr2(find(abs(perr2)>200))=nan;

figure(3);clf

h1=histogram(perr1(:),30,'Normalization','pdf');
hold on
h2=histogram(perr2(:),h1.BinEdges,'Normalization','pdf');
freqline(nanmean(perr1),'b--')
freqline(nanmean(perr2),'r--')
xlabel('% Error in \chi','fontsize',16)
ylabel('Pdf ','fontsize',16)
title('Error From Neglect of Horizontal Velocity')
grid on
xlim(100*[-1 1])
set(gca,'Xtick',[-100:25:100])
legend([h1 h2],'+0.1m/s','+1m/s')
% vline(nanmedian(perr1),'b--')
% vline(nanmedian(perr2),'r--')
%legend([h2],'+0.3m/s')

text(25,1.25e4,['\mu=' num2str(roundx(nanmean(perr1),2))],'color','b','fontsize',16)
text(25,1e4,['\mu=' num2str(roundx(nanmean(perr2),2))],'color','r','fontsize',16)

text(60,1.25e4,['\sigma=' num2str(roundx(nanstd(perr1),2))],'color','b','fontsize',16)
text(60,1e4,['\sigma=' num2str(roundx(nanstd(perr2),2))],'color','r','fontsize',16)
%%
if saveplot==1
    SetPaperFigPath
print(fullfile(figdir,'Hist_perr_fspdvary'),'-dpng')
end
%%