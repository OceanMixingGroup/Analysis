%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotExampleSpec.m
%
% * Makes plot for chipod methods paper *
%
%
%---------------
% 04/11/16 - A.Pickering - apickering@coas.oregonstate.edu
% 04/15/16 - AP - Remove cham chi where cham eps=NaN, add eps values
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
%
% clear ; close all
%
% datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'
%
% saveplots=1
% castnum=2700
%
% fmax=7
%
% load( fullfile(datdir,['zsm10m_fmax' num2str(fmax) 'Hz_respcorr0_fc_99hz_gamma20'],['EQ14_' sprintf('%04d',castnum) 'avg.mat']) )
% avg1=avg;clear avg
%
% % load chameleon profile also to compare
%
% load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',castnum) '.mat'])
% cham=avg;clear avg
%
% % don't use chis where epsilon is NaN
% idb=find(isnan(cham.EPSILON));
% cham.CHI(idb)=nan;
%
% chameps=interp1(cham.P,cham.EPSILON,avg1.P);
% chamchi=interp1(cham.P,cham.CHI,avg1.P);
%
% %%
%
% for iwind=1:5:length(avg1.P)
%
%     ax=PlotExSpec(avg1,chamchi,chameps,iwind)
%
%     %     figure(1);clf
%     %
%     %     h1=loglog(avg1.ks(iwind,:),avg1.kspec(iwind,:),'k');
%     %     hold on
%     %     loglog(avg1.kks(iwind,:),avg1.kkspec(iwind,:),'k--');
%     %
%     %     % compute kraichnan spec for actual value measured by chameleon also
%     %     clear nu tdif k chi eps kb spec_vals
%     %     nu=avg1.nu(iwind);
%     %     tdif=avg1.tdif(iwind);
%     %     qq=7;
%     %     k=avg1.kks(iwind,:);
%     %     chi=chamchi(iwind);
%     %     eps=chameps(iwind);
%     %     kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
%     %     [spec_vals]=kraichnan(nu,k,kb,tdif,chi,qq);
%     %     hcham=loglog(k,spec_vals,'m','linewidth',2)
%     %     %
%     %     ylim([1e-12 1e0])
%     %     grid on
%     %     xlabel('wavenumber [cpm]','fontsize',16)
%     %     ylabel('\Phi_{dT/dz}[K^2m^{-2}cpm^{-1}]','fontsize',16)
%     %     xlim([0.1 1.5e3])
%     %
%     %     % mark kmax
%     %     freqline(avg1.Params.fmax./avg1.fspd(iwind),'r--')
%     %     text(avg1.Params.fmax./avg1.fspd(iwind),1e-2,'k_{max}')
%     %
%     %     % mark batchelor wavenumber k_b
%     %     freqline(kb,'c--')
%     %     text(kb,1e-2,'k_b')
%     %
%     %     title(['P=' num2str(avg1.P(iwind)) 'db, chameleon \epsilon=' num2str(chameps(iwind)) ])
%     %     %legend([h1 h2 h3],['fc=' num2str(avg1.Params.fc)],['fc=' num2str(avg2.Params.fc)],['fc=' num2str(avg3.Params.fc)])
%     %     % legend([h1 h2 hcham],'\chi pod corr.','\chi pod no corr.','chameleon')
%     %     legend([h1 hcham],'\chi pod ','chameleon','location','northwest')
%     %
%     %     freqline(avg1.kstart(iwind),'b--')
%     %     freqline(avg1.kstop(iwind),'b--')
%     %
%     %     text(1e-1,1e-8,['\chi_{\chi} =' sprintf('%.2e',avg1.chi1(iwind))],'fontweight','bold')
%     %     text(1e-1,1e-9,['\chi_{\epsilon} =' sprintf('%.2e',chamchi(iwind))],'color','m','fontweight','bold')
%     %     text(1e-1,1e-10,['\chi_{\chi}/\chi_{\epsilon} =' sprintf('%.2f',avg1.chi1(iwind)/chamchi(iwind))],'color','b','fontweight','bold')
%     %
%     %     text(1e0,1e-8,['\epsilon_{\chi} =' sprintf('%.2e',avg1.eps1(iwind))],'fontweight','bold')
%     %     text(1e0,1e-9,['\epsilon_{\epsilon} =' sprintf('%.2e',chameps(iwind))],'color','m','fontweight','bold')
%     %     text(1e0,1e-10,['\epsilon_{\chi}/\epsilon_{\epsilon} =' sprintf('%.2f',avg1.eps1(iwind)/chameps(iwind))],'color','b','fontweight','bold')
%
%     pause(0.2)
%
%     if saveplots==1
%         %        figdir='/Users/Andy/Dropbox/AP_Share_With_JN/chipodspec'
%         figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Figures/ExSpec/'
%         print( fullfile(figdir,['wind' num2str(iwind) '_specs']),'-dpng')
%     end
%
% end



%% Look at spectra for only certain epsilon ranges

%[val,I]=nanmax(chameps)
%[val,I]=nanmax(chamchi)
%ig=find(log10(chameps)>-6)
%ig=find(log10(chamchi)>-7)
%ig=find(log10(chamchi)<-9)
%ig=isin(log10(chamchi),[-10 -9])
%ig=isin(log10(chamchi),[-9 -8])
%ig=isin(log10(chamchi),[-8 -7])
%ig=isin(log10(chamchi),[-7 -6])
%ig=find(log10(chameps)>-6)

%
% for whi=1:length(ig)
%
%     whi
%     I=ig(whi)
%
%     iwind=I
%     figure(1);clf
%     ax=PlotExSpec(avg1,chamchi,chameps,iwind)
%
%     pause%(0.1)
%
% end

%%

% ** Makes plot currently used in paper - June 3 2016 **

clear ; close all

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'

saveplots=1
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
ax1=PlotExSpec(avg1,chamchi,chameps,iwind)
%pause
%end
%
subplot(212)
whi=5
I=ig2(whi)
iwind=I
ax2=PlotExSpec(avg1,chamchi,chameps,iwind)

%linkaxes([ax1 ax2])

if saveplots==1
    SetPaperFigPath
    print(fullfile(figdir,'ExSpec_HiLowEps'),'-dpng')
end
%%
