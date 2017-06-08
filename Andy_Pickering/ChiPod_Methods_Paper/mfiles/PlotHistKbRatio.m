%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotHistKbRatio.m
%
% * Makes plot for chipod methods paper *
%
% Make a histogram of the ratio of kmax to kb for the EQ14 data, showing
% how much of the spectrum we resolve.
%
%
%---------------
% 04/15/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'

saveplots=1

kbrat=[]
kball=[]
kmaxall=[]
fspdall=[]
fmax=7
hb=waitbar(0)
for castnum=1:3100
    waitbar(castnum/3100,hb)
    %castnum=7
    clear avg1 avg cham
    
    try
        
        load( fullfile(datdir,['zsm10m_fmax' num2str(fmax) 'Hz_respcorr0_fc_99hz_gamma20'],['EQ14_' sprintf('%04d',castnum) 'avg.mat']) )
        avg1=avg;clear avg
        
        
        % load chameleon profile
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',castnum) '.mat'])
        cham=avg;clear avg
        
        % don't use chis where epsilon is NaN
        clear idb chameps
        idb=find(isnan(cham.EPSILON));
        idb=find(log10(cham.EPSILON)>-5);
        cham.EPSILON(idb)=nan;
        idb=find(log10(cham.EPSILON)<-12);
        cham.EPSILON(idb)=nan;
        %cham.CHI(idb)=nan;
        
        chameps=interp1(cham.P,cham.EPSILON,avg1.P);
        
        clear nu tdif qq eps kb
        nu=avg1.nu;%
        tdif=avg1.tdif;%
        qq=7;
        eps=chameps;%
        kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
        %
        clear kmax fspd
        fspd=avg1.fspd;
        fspd( find( abs(fspd)<0.1) )=nan;
        
        kmax=fmax./abs(fspd);
        kmaxall=[kmaxall ; kmax];
        kball=[kball ; kb];
        kbrat=[kbrat ; kmax./kb];
        fspdall=[fspdall ; fspd];
        
    end
end
delete(hb)

%%

ib=find(abs(fspdall)>1.5);
kbrat(ib)=nan;
fspdall(ib)=nan;

ib=find(abs(fspdall)<0.2);
kbrat(ib)=nan;
fspdall(ib)=nan;

Nm='pdf'

figure(1);clf

subplot(121)
histogram(kbrat(:),'edgecolor','none','Normalization',Nm)
xlim([0 0.6])
xlabel('k_{max}/k_b','fontsize',16)
ylabel(Nm,'fontsize',16)
grid on
title('Ratio of k_{max} / k_b ')

subplot(122)
histogram(abs(fspdall(:)),'edgecolor','none','Normalization',Nm)
xlabel('u [ms^{-1}]','fontsize',16)
ylabel(Nm,'fontsize',16)
grid on

if saveplots==1
    SetPaperFigPath
print(fullfile(figdir,'Hist_kbrat'),'-dpng')
end

%%