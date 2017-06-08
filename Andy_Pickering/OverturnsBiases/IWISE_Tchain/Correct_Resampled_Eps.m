%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Correct_Resampled_Eps.m
%
% Try to identify and correct biased epsilon in resampled T-chain data by
% examing turbulent velocity scale.
%
% Modified from part of Misc16Mar.m
%----------------
% 23 November 2015 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Use Lt*N, which is a vertical velocity scale for the overturn?

clear ; close all

saveplots=1

whmoor=3; minOT=50; 
%testnum=3 % 2=0.15
%testnum=2 % w=0.5
testnum=1 % w=0.25
%testnum=4 % w=0.75


cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true data'
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
%


wt_true=xx2.Lt_each.*sqrt(xx2.Otnsq_each);

eps1=[]
eps2=[]
wt2=[];
for whc=1:REsamp.Nshift
    clear wt_samp ib
%   wt_samp=REsamp.Lot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
    wt_samp=REsamp.Lttot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
    wt_samp=real(wt_samp);
    wt2=[wt2 ; wt_samp(:)];
    
    clear ec
    ec=squeeze(REsamp.eps(:,:,whc));
    eps1=[eps1 nanmean(ec,2)];
    ib=find(wt_samp> (0.5*REsamp.w_samp)  );
    ec(ib)=nan;
    eps2=[eps2 nanmean(ec,2)];
    
end
%
Nm='pdf'
Ds='stair'
figure(1);clf
ht=histogram(wt_true,'Normalization',Nm,'DisplayStyle',Ds)
hold on
%hs=histogram(wt_samp,'BinEdges',ht.BinEdges,'Normalization',Nm,'DisplayStyle',Ds)
hs=histogram(wt2,'BinEdges',ht.BinEdges,'Normalization',Nm,'DisplayStyle',Ds)
legend('true','samp')
freqline(REsamp.w_samp)
ylabel(Nm)
xlabel('w_{turb}')
grid on

%% Nan out resampled epsilon where w_t>w_samp, and plot mean profile

figure(2);clf
agutwocolumn(0.7)
wysiwyg
semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)
hold on
semilogx(nanmean(eps1,2),REsamp.z,'linewidth',2)
semilogx(nanmean(eps2,2),REsamp.z,'linewidth',2)
xlim([1e-9 1e-6])
ylim([600 nanmax(xx2.z)])
axis ij 
grid on
legend('true','samp','samp corr')
xlabel('<\epsilon> [Wkg^{-1}]','fontsize',16)
ylabel('Depth [m]','fontsize',16)
title(['Tchain ' num2str(whmoor) ' - W_{samp}=' num2str(REsamp.w_samp)])

if saveplots==1
   figname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_Testnum' num2str(testnum) '_Raw_corr_profiles']) 
   print(figname,'-dpng')
end
%%