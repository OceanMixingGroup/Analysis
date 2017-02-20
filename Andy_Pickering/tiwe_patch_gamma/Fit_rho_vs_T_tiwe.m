%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Fit_rho_vs_T.m
%
%
% Fit a line to rho vs T to get constant of proportionality alpha
% ie sgth=alpha*theta
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% use binned data?
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/TIWE/data/tiwe_1mavg_combined.mat')
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/

cham.THETA=sw_ptmp(cham.S,cham.T,cham.P,0);
cham.SIGMA=sw_pden(cham.S,cham.T,cham.P,0);

x=cham.THETA(:);
y=cham.SIGMA(:);

ig=find(~isnan(x) & ~isnan(y) );

P=polyfit(x(ig),y(ig),1);

figure(1);clf
plot(x(ig),y(ig),'.','color',0.75*[1 1 1])
ylim([1022 1027])
grid on
hold on
plot(x(ig),polyval(P,x(ig)),'linewidth',2)
xlabel('theta')
ylabel('sgth')

%%
SetNotesFigDir
print( fullfile(NotesFigDir,'tiwe_sght_vs_theta_fit'), '-dpng')

%% Do we get same answer using in-situ density/temp?

x=cham.T(:);
y=cham.SIGMA(:);

ig=find(~isnan(x) & ~isnan(y) );

P=polyfit(x(ig),y(ig),1);

figure(1);clf
plot(x(ig),y(ig),'.')
ylim([1022 1027])
grid on
hold on
plot(x(ig),polyval(P,x(ig)),'linewidth',2)
xlabel('theta')
ylabel('sgth')

P(1)


%% Does answer change if we use raw (isntead of 1-m avg) data?
% 
% % no, alpha is nearly the same..
% 
% clear ; close all
% 
% datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'
% 
% theta=[];
% sgth=[];
% 
% for cnum=1:50:3100
%     
%     disp(['Working on cnum ' num2str(cnum) ])
%     
%     close all
%     clear cal head
%     clear tpspec kspec kkspec fspec kks ks
%     
%     try
%         %%
%         % Load the data for this cast
%         load( fullfile(datdir,['EQ14_' sprintf('%04d',cnum) '.mat']) )
%         
%         clear cal
%         cal=cal2;
%         clear cal2
%         
%         
%         theta=[theta(:) ; sw_ptmp(cal.SAL,cal.T1,cal.P,0)];
%         sgth=[sgth(:) ; sw_pden(cal.SAL,cal.T1,cal.P,0)];
%         
%     end % try
%     
%     
% end % cnum
% %
% 
% figure(1);clf
% loglog(theta,sgth+1001,'.')
% grid on
% 
% x=theta;
% y=sgth+1000;
% ig=find(~isnan(x) & ~isnan(y) );
% 
% P=polyfit(x(ig),y(ig),1);
% P(1)
% %%
