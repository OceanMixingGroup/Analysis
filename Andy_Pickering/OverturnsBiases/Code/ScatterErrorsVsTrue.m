%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ScatterErrorsVsTrue.m
%
% Modified from parts of PlotResampVsTrue_Direct.m
%
%--------------
% 01/25/15 - A.Pickering
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplot=0
whmoor=3; xl=500*[-1 1]
testnum=1
minOT=50

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%~~ load 1 resampled dataset to get sizes etc.
whcase=1
clear fname x_resamp
%fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
load(fullfile('Data/',fname))

%~~  Load full T-chain data
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

%~~ load each resamp case (shifted in time) and save true data
% corresponding to those times

T_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
T_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

clear Lot_Allcases_Resamp Lot_Allcases_True
Lot_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
Lot_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

% add Lttot (rms thorpe scales)

Lttot_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
Lttot_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case


Eps_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
Eps_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

for whcase=1:100
    clear fname x_resamp
    %    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
    load(fullfile('Data/',fname))
    
    % find indices of Tchain data corresponding to time of each profile
    clear It
    It=nan*ones(size(x_resamp.time));
    for wht=1:length(x_resamp.time)
        [val,It(wht)]=nanmin(abs(xx2.yday-x_resamp.time(wht))) ;
        if val>1/24
            It(wht)=nan;
        end
    end
    
    %~
    
    Lot_Allcases_Resamp(:,whcase)=x_resamp.Lot(:);
    Lttot_Allcases_Resamp(:,whcase)=x_resamp.Lttot(:);
    
    clear Ltemp
    Ltemp=xx2.Lot(:,It);
    Lot_Allcases_True(:,whcase)=Ltemp(:);
    
    clear Ltemp
    Ltemp=xx2.Lttot(:,It);
    Lttot_Allcases_True(:,whcase)=Ltemp(:);
    
    clear Epstemp
    Eps_Allcases_Resamp(:,whcase)=x_resamp.eps(:);
    
    Epstemp=xx2.eps(:,It);
    Eps_Allcases_True(:,whcase)=Epstemp(:);
    
    %`
    clear Ttemp
    T_Allcases_Resamp(:,whcase)=x_resamp.T(:);
    
    Ttemp=xx2.T(:,It);
    T_Allcases_True(:,whcase)=Ttemp(:);
    
    
end % whcase

% Scatter plot the difference in Lt vs actual Lt

clear LtDif X Y ig P x y
% use individual overturn sizes
LtDif=Lot_Allcases_Resamp(:)-Lot_Allcases_True(:);

% use rms thorpe scale
%LtDif=Lttot_Allcases_Resamp(:)-Lttot_Allcases_True(:);

figure(1);clf
agutwocolumn(1)
ax = MySubplot(0.12, 0.03, 0.13, 0.06, 0.1, 0.1, 2,2);
wysiwyg

axes(ax(1))
%plot(Lot_Allcases_True(:),abs(LtDif),'o')
%plot(Lot_Allcases_True(:),(LtDif),'.')
plot(Lttot_Allcases_True(:),(LtDif),'.')
grid on
xlabel('L_{true} (m) ','fontsize',16)
ylabel('L_{resamp}-L_{true} (m) ','fontsize',16)

clear X Y ig P x y
X=Lttot_Allcases_True(:);
Y=LtDif(:);
ig=find(~isnan(X) & ~isnan(Y));
X=X(ig);Y=Y(ig);
P=polyfit(X,Y,1)
x=linspace(nanmin(X),nanmax(X),100)
y=polyval(P,x);
hold on
plot(x,y,'.-')
xlim([0 350])
ylim(600*[-1 1])
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])

%
% Scatter plot the difference in Eps vs actual Eps

clear LtDif X Y ig P x y
% use individual overturn sizes
EpsDif=Eps_Allcases_Resamp(:)-Eps_Allcases_True(:);

EpsDif(log10(Eps_Allcases_True)<-10)=nan;
%
ip=find(EpsDif>0);
in=find(EpsDif<0);

axes(ax(2))
hp=loglog(Eps_Allcases_True(ip),EpsDif(ip),'.')
hold on
hn=loglog(Eps_Allcases_True(in),abs(EpsDif(in)),'.','color',0.7*[1 1 1])
xlim([1e-10 1e-3])
ylim([1e-10 1e-3])
grid on
xlabel('\epsilon _{true} (Wkg^{-1}) ','fontsize',16)
ylabel('\epsilon _{resamp}-\epsilon _{true} (Wkg^{-1}) ','fontsize',16)
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])
legend([hp hn],'bias>0','bias<0','location','best')

% Scatter plot the difference in Eps vs Lt

clear LtDif X Y ig P x y
% use individual overturn sizes
EpsDif=Eps_Allcases_Resamp(:)-Eps_Allcases_True(:);

EpsDif(log10(Eps_Allcases_True)<-10)=nan;
%
ip=find(EpsDif>0);
in=find(EpsDif<0);

axes(ax(3))
hp=loglog(Lot_Allcases_True(ip),EpsDif(ip),'.')
hold on
hn=loglog(Lot_Allcases_True(in),abs(EpsDif(in)),'.','color',0.7*[1 1 1])
xlim([10^(1.4) 1e3])
ylim([1e-10 1e-3])
xlabel('L_{true} (m) ','fontsize',16)
ylabel('\epsilon _{resamp}-\epsilon _{true} (Wkg^{-1}) ','fontsize',16)
grid on ; gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])
legend([hp hn],'bias>0','bias<0','location','best')
%
clear X Y ig P x y
X=Lot_Allcases_True(ip);
Y=abs(EpsDif(ip));
ig=find(~isnan(X) & ~isnan(Y));
X=X(ig);Y=Y(ig);
P=polyfit(X,Y,1)
x=linspace(nanmin(X),nanmax(X),100);
y=polyval(P,x);
hold on
%
loglog(x,y,'.-','color','k')
shg
%
% clear X Y ig P x y
% X=Lot_Allcases_True(:);
% Y=abs(EpsDif(:));
% ig=find(~isnan(X) & ~isnan(Y));
% X=X(ig);Y=Y(ig);
% P=polyfit(X,Y,1)
% x=linspace(nanmin(X),nanmax(X),100);
% y=polyval(P,x);
% hold on
% %
% loglog(x,y,'.-','color','r')
% shg
% %
% clear X Y ig P x y
% X=Lot_Allcases_True(in);
% Y=abs(EpsDif(in));
% ig=find(~isnan(X) & ~isnan(Y));
% X=X(ig);Y=Y(ig);
% P=polyfit(X,Y,1)
% x=linspace(nanmin(X),nanmax(X),300);
% y=polyval(P,x);
% hold on
% %
% loglog(x,y,'.-','color','m')
%shg
%loglog(X,0.64*(X.^2)./2e13,'.-','color','y')
%

TDif=T_Allcases_Resamp(:)-T_Allcases_True(:);

axes(ax(4))
plot(T_Allcases_True(:),TDif(:),'.')
grid on
xlabel('T_{true} (m) ','fontsize',16)
ylabel('T_{resamp}-T_{true} (m) ','fontsize',16)

% X=T_Allcases_True(:);
% Y=TDif(:);
% ig=find(~isnan(X) & ~isnan(Y));
% X=X(ig);Y=Y(ig);
% P=polyfit(X,Y,1)
% x=linspace(nanmin(X),nanmax(X),100)
% y=polyval(P,x);
% hold on
% plot(x,y,'.-')
xlim([2 8])
ylim(2*[-1 1])
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])
%
%%
if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_Testnum' num2str(testnum) '_ScatterErrorVsTrue'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf','-r300')
end

%%