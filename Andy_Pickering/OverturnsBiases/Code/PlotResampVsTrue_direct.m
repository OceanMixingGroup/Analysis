%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotResampVsTrue_direct.m
%
% Plot resampled versus true overturn sizes and epsilon
% from T-chains, for data at same time and depths.
%
% Resampling done in ResampMPnew12Nov.m
%
% Use 'true' proifles at times corresponding to midpoint of resampled
% profiles.
%
% See also PlotResampResults12Nov.m
%
%
% Makes 4 panel plot, plot hist of LtDif for different speeds on
% same x-axes. I think all will be gaussian but get narrower with
% increasing speed.
%
%-------------
% 12/17/14  - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Mmake 4 panel plot, plot hist of LtDif for different speeds on
% same x-axes. I think all will be gaussian but get narrower with
% increasing speed.
%

clear ; close all

saveplot=0
whmoor=3; xl=500*[-1 1]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.09, 2,2);

whax=1;

for testnum=4%[3 1 2 4]
    
    %~~ load 1 resampled dataset to get sizes etc.
    whcase=1
    clear fname x_resamp
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
    load(fullfile('Data/',fname))
    
    %~~  Load full T-chain data
    load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))
    
    %~~ load each resamp case (shifted in time) and save true data
    % corresponding to those times
    
    % T_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
    % T_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
    
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
        fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
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
        %     clear Ttemp
        %     T_Allcases_Resamp(:,whcase)=x_resamp.T(:);
        
        %     Ttemp=xx2.T(:,It);
        %     T_Allcases_True(:,whcase)=Ttemp(:);
        
        
    end % whcase
    
    clear LtDif
    % use individual overturn sizes
    LtDif=Lot_Allcases_Resamp(:)-Lot_Allcases_True(:);
    % use rms thorpe scale
    %    LtDif=Lttot_Allcases_Resamp(:)-Lttot_Allcases_True(:);
    
    axes(ax(whax))
    h=histogram(LtDif,[-600:25:600])
    xlim(xl)
    
    % fit a normal distribution
    XGrid=linspace(-500,500,100);
    pd1 = fitdist(LtDif(:), 'normal');
    YPlot = pdf(pd1,XGrid);
    lims=axis;
    ymax=lims(4);
    scale=ymax/max(YPlot);
    hold on
    hLine = plot(XGrid,YPlot*scale,'r','LineWidth',2)
    grid on
    ylabel('#')
    xlabel('L_{resamp} - L_{true}')
    %    clear EpsDif
    %   EpsDif=Eps_Allcases_Resamp(:)-Eps_Allcases_True(:);
    title(['\mu =' num2str(roundx(nanmean(LtDif),2)) ', \sigma=' num2str(roundx(nanstd(LtDif),2))])
    SubplotLetterMW(['w=' num2str(x_resamp.w_samp) 'ms^{-1}'],0.03,0.05)
    shg
    whax=whax+1
    
end % testnum

linkaxes(ax,'x')

if saveplot==1
    %
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_LtDifHist_4Cases'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
    %
end

%%

Nm='probability'
Nm='count'
figure(33);clf
histogram(Lot_Allcases_Resamp,'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(Lot_Allcases_True,'DisplayStyle','Stair','Normalization',Nm)
%%
%
%
%% Scatter plot the difference in Lt vs actual Lt

clear LtDif X Y x y ig P
%LtDif=Lot_Allcases_Resamp(:)-Lot_Allcases_True(:);
LtDif=Lttot_Allcases_Resamp(:)-Lttot_Allcases_True(:);

figure(1);clf
agutwocolumn(0.6)
wysiwyg
%plot(Lot_Allcases_True(:),abs(LtDif),'o')
%plot(Lot_Allcases_True(:),(LtDif),'.')
plot(Lttot_Allcases_True(:),(LtDif),'.')
grid on
xlabel('L_{true} (m) ','fontsize',16)
ylabel('L_{resamp}-L_{true} (m) ','fontsize',16)

X=Lttot_Allcases_True(:);
Y=LtDif(:);
ig=find(~isnan(X) & ~isnan(Y));
X=X(ig);Y=Y(ig);
P=polyfit(X,Y,1)
x=linspace(nanmin(X),nanmax(X),100)
y=polyval(P,x);
hold on
plot(x,y,'.-')
%xlim([0 550])
%ylim(275*[-1 1])
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])

%%
clear EpsDif X Y x y ig P
EpsDif=Eps_Allcases_Resamp(:)-Eps_Allcases_True(:);

ib=find(log10(Eps_Allcases_True(:))<-10);
EpsDif(ib)=nan;

figure(2);clf
agutwocolumn(0.6)
wysiwyg
%plot(Lot_Allcases_True(:),abs(LtDif),'o')
loglog(Eps_Allcases_True(:),(EpsDif),'.')
grid on
xlabel('\epsilon_{true}','fontsize',16)
ylabel('\epsilon_{resamp}-\epsilon_{true}','fontsize',16)
%
% X=Eps_Allcases_True(:);
% Y=EpsDif(:);
% ig=find(~isnan(X) & ~isnan(Y));
% X=X(ig);Y=Y(ig);
% P=polyfit(X,Y,1)
% x=linspace(nanmin(X),nanmax(X),100)
% y=polyval(P,x);
% hold on
% plot(x,y,'.-')
%xlim([0 550])
%ylim(275*[-1 1])
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])
%


clear EpsDif X Y x y ig P
EpsDif=Eps_Allcases_Resamp(:)-Eps_Allcases_True(:);

%ib=find(log10(Eps_Allcases_True(:))<-10);
%EpsDif(ib)=nan;

figure(3);clf
agutwocolumn(0.6)
wysiwyg
%plot(Lot_Allcases_True(:),abs(LtDif),'o')
plot(Lot_Allcases_True(:),(EpsDif),'.')
%plot(Lttot_Allcases_True(:),(EpsDif),'.')
grid on
xlabel('L_{true}','fontsize',16)
ylabel('\epsilon_{resamp}-\epsilon_{true}','fontsize',16)

X=Lttot_Allcases_True(:);
Y=EpsDif(:);
ig=find(~isnan(X) & ~isnan(Y));
X=X(ig);Y=Y(ig);
P=polyfit(X,Y,1)
x=linspace(nanmin(X),nanmax(X),100)
y=polyval(P,x);
hold on
plot(x,y,'.-')
%xlim([0 550])
%ylim(275*[-1 1])
gridxy
title(['T-chain ' num2str(whmoor) ', w=' num2str(x_resamp.w_samp) 'm/s'])

%%

clear ; %close all

saveplot=0

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3;
testnum=4

%~~ load 1 resampled dataset to get sizes etc.
whcase=1
clear fname x_resamp
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
load(fullfile('Data/',fname))

%~~  Load full T-chain data
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps']))

%~~ load each resamp case (shifted in time) and save true data
%correspondign to those times

T_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
T_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

Lot_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
Lot_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

Eps_Allcases_Resamp=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case
Eps_Allcases_True=nan*ones(length(x_resamp.time)*length(x_resamp.z),100); % each profile is one case

for whcase=1:100
    clear fname x_resamp
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
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
    clear Lot_all_resamp Ltemp
    Lot_Allcases_Resamp(:,whcase)=x_resamp.Lot(:);
    
    Ltemp=xx2.Lot(:,It);
    Lot_Allcases_True(:,whcase)=Ltemp(:);
    
    %`
    clear Ttemp
    T_Allcases_Resamp(:,whcase)=x_resamp.T(:);
    
    Ttemp=xx2.T(:,It);
    T_Allcases_True(:,whcase)=Ttemp(:);
    
    %~
    clear Epstemp
    Eps_Allcases_Resamp(:,whcase)=x_resamp.eps(:);
    
    Epstemp=xx2.eps(:,It);
    Eps_Allcases_True(:,whcase)=Epstemp(:);
    
    
    %     figure(1);clf
    %     plot(Ltemp(:),x_resamp.Lot(:),'o')
    %     title(['case ' num2str(whcase)])
    %     pause(0.1)
end

%

LtDif=Lot_Allcases_Resamp(:)-Lot_Allcases_True(:);

%
figure;clf
agutwocolumn(0.5)
wysiwyg
%axes(ax(4))
h=histogram(LtDif,50)
xlim(nanmax(abs(LtDif))*[-1 1])
grid on
ylabel('#')
xlabel('L_{resamp} - L_{true}')
SubplotLetterMW(['Tchain ' num2str(whmoor)])
shg
%

if saveplot==1
    %
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_LtDifHist'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
    %
end
%%

saveplot=0

figure(1);clf

agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.09, 2,2);

axes(ax(1))
%loglog(Lot_Allcases_True(:),Lot_Allcases_Resamp(:),'o')
plot(Lot_Allcases_True(:),Lot_Allcases_Resamp(:),'.')
%axis equal
xvec=0:10:nanmax(Lot_Allcases_True(:));
hold on
plot(xvec,xvec,'k')
grid on
SubplotLetterMW('L_{OT}')
xlabel('true')
ylabel('resamp')
%

%figure(1);clf
axes(ax(2))
%loglog(Lot_Allcases_True(:),Lot_Allcases_Resamp(:),'o')
plot(T_Allcases_True(:),T_Allcases_Resamp(:),'.')
%axis equal
xvec=2:8;
hold on
plot(xvec,xvec,'k')
grid on
SubplotLetterMW('T')
xlabel('true')
ylabel('resamp')

%

%figure(1);clf
axes(ax(3))
plot(log10(Eps_Allcases_True(:)),log10(Eps_Allcases_Resamp(:)),'.')
%axis equal
xvec=-10:-3;
hold on
plot(xvec,xvec,'k')
grid on
%ylabel('log_{10} \epsilon Resamp')
%xlabel('log_{10} \epsilon True')
SubplotLetterMW('log_{10} \epsilon')
xlabel('true')
ylabel('resamp')

log10(nanmean(Eps_Allcases_Resamp(:)))
log10(nanmean(Eps_Allcases_True(:)))
%

%figure(4);clf
axes(ax(4))
hist(Lot_Allcases_Resamp(:)-Lot_Allcases_True(:),50)
xlim(600*[-1 1])
grid on
ylabel('#')
xlabel('L_{resamp} - L_{true}')



if saveplot==1
    %
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_DirectCompare'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
    %
end

%%
figure(4);clf
hist(T_Allcases_Resamp(:)-T_Allcases_True(:),100)
xlim(0.5*[-1 1])
grid on
ylabel('#')
xlabel('T_{resamp} - T_{true}')

%%
figure(4);clf
hist( log10(Eps_Allcases_Resamp(:)) - log10(Eps_Allcases_True(:)) ,10)
%xlim(600*[-1 1])
grid on
ylabel('#')
%xlabel('L_{resamp} - L_{true}')

%%

figure(1);clf
plot(x_resamp.time,xx2.yday(It),'o-')
grid on

T_all_resamp=x_resamp.T(:);
Ttemp=xx2.T(:,It);
T_all_true=Ttemp(:);

Lot_all_resamp=x_resamp.Lot(:);
Ltemp=xx2.Lot(:,It);
Lot_all_true=Ltemp(:);

eps_all_resamp=x_resamp.eps(:);
etemp=xx2.eps(:,It);
eps_all_true=etemp(:);

figure(1);clf
%plot(Lot_all_true,Lot_all_resamp,'o')
loglog(T_all_true,T_all_resamp,'d')
grid on
xvec=linspace(nanmin(T_all_true),nanmax(T_all_true),100);
hold on
plot(xvec,xvec,'k')


figure(2);clf
%plot(Lot_all_true,Lot_all_resamp,'o')
loglog(Lot_all_true,Lot_all_resamp,'d')
grid on
xvec=linspace(nanmin(Lot_all_true),nanmax(Lot_all_true),100);
hold on
plot(xvec,xvec,'k')


figure(3);clf
%plot(Lot_all_true,Lot_all_resamp,'o')
loglog(eps_all_true,eps_all_resamp,'d')
grid on

xvec=linspace(nanmin(eps_all_true),nanmax(eps_all_true),100);
hold on
plot(xvec,xvec,'k')
%%