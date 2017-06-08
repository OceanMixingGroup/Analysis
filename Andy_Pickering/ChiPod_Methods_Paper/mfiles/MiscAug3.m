%%
% ** Makes plot currently used in paper - June 3 2016 **

clear ; close all

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'

%saveplots=1
castnum=2700

fmax=7

load( fullfile(datdir,['zsm10m_fmax' num2str(fmax) 'Hz_respcorr0_fc_99hz_gamma20'],['EQ14_' sprintf('%04d',castnum) 'avg.mat']) )
avg1=avg;clear avg

% load chameleon profile also to compare

%load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',castnum) '.mat'])
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',castnum) '.mat'])
cham=avg;clear avg

% don't use chis where epsilon is NaN
idb=find(isnan(cham.EPSILON));
cham.CHI(idb)=nan;

chameps=interp1(cham.P,cham.EPSILON,avg1.P);
chamchi=interp1(cham.P,cham.CHI,avg1.P);
%%

ig1=find(log10(chameps)<-7.5 & log10(chameps)>-8.5 & (avg1.eps1./ chameps) >0.1 )
ig2=find(log10(chamchi)>-6 )

figure(1);clf
agutwocolumn(1)
wysiwyg
subplot(211)

whi=25
%for whi=24:length(ig1)
I=ig1(whi)
iwind=I
ax=PlotExSpec(avg1,chamchi,chameps,iwind)
%pause
%end
%
subplot(212)
whi=7
I=ig2(whi)
iwind=I
ax=PlotExSpec(avg1,chamchi,chameps,iwind)

% 
% if saveplots==1
%     SetPaperFigPath
%     print(fullfile(figdir,'ExSpec_HiLowEps'),'-dpng')
% end
%$
ig2=find(log10(chameps)>-7)
for whi=1:length(ig2)
%    whi=7
I=ig2(whi)
iwind=I
figure(1);clf
agutwocolumn(0.6)
wysiwyg
ax=PlotExSpec(avg1,chamchi,chameps,iwind)
pause
end

%%

ig1=find(log10(chameps)<-8.5)
k_1=avg1.ks(ig1,:);
ki_1=k_1(1,:);
kspec1=avg1.kspec(ig1,:);

kspeci_1=nan*ones(length(ig1),length(ki_1));
for whz=1:length(ig1)
kspeci_1(whz,:)=interp1(k_1(whz,:),kspec1(whz,:),ki_1);   
end


ig2=find(log10(chameps)>-6.5)
k_2=avg1.ks(ig2,:);
ki_2=k_2(1,:);
kspec2=avg1.kspec(ig2,:);

kspeci_2=nan*ones(length(ig2),length(ki_2));
for whz=1:length(ig2)
kspeci_2(whz,:)=interp1(k_2(whz,:),kspec2(whz,:),ki_2);   
end
%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
loglog(ki_1,(kspeci_1))
hold on
loglog(ki_1,nanmean(kspeci_1),'k','linewidth',3)
%loglog(ki_2,nanmean(kspeci_2))
grid on
ylim([1e-13 1e0])
xlim([1e0 1e2])

subplot(212)
loglog(ki_2,(kspeci_2))
hold on
loglog(ki_2,nanmean(kspeci_2),'k','linewidth',3)
%loglog(ki_2,nanmean(kspeci_2))
grid on
ylim([1e-13 1e0])
xlim([1e0 1e2])
%%

figure(2);clf

loglog(ki_1,nanmean(kspeci_1))
hold on
loglog(ki_2,nanmean(kspeci_2))
grid on
legend('1','2')
freqline(2)
freqline(7)
%loglog(avg1.ks(ig2,:),avg1.kspec(ig2,:),'.')
%hold on
%loglog(ki,kspeci)
%hold on
%loglog(ki,nanmean(kspeci),'k','linewidth',2)

%%