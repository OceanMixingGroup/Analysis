%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% EvalFalseOverturns.m
%
% Examine effect of false overturns (overturns in resampled data where real
% profile was stable) on depth and time averages of epsilon etc.
%
%
% 28 Jan 2015 - A. Pickering
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Contour eps and isotherms, and plot where overturns are introduced into
% profiles that were actually stable.
% 27 Jan 2015
% 28 Jan 2015

clear ; close all

plotiso=1
plotsamp=1
whmoor=3 % mooring #
testnum=1
minOT=50

xl=[168 170]
xl=[168.3 168.4]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
%xl=[183.2 183.7]
%xl=[210 212]
%xl=[165 185]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true' T-chain data
load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

figure(2);clf
agutwocolumn(0.6)
wysiwyg
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));
id=isin(xx2.yday,xl);

ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14
%
caxis(cl)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ],'interpreter','none')

if plotiso==1
    hold on
    %contour(xx2.yday(id),xx2.z,xx2.T(:,id),[1:7],'k')
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:10:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)

%  load resampled data
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%
LotAll_real=[];
LotAll_samp=[];

for whc=80%:100
    
    clear idd tvec
    
    % times of midpoints of resampled profiles
    tvec=REsamp.timeall(whc,:);
    
    % find profiles in the time range we are looking at
    idd=isin(tvec,xl);
    
    % plot resampling path (depth-time)
    if plotsamp==1
        plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')
    end
    
    clear Itt
    Itt=nan*ones(length(idd),1);
    for tt=1:length(idd)
        [val,Itt(tt)]=nanmin( abs( tvec(idd(tt)) - xx2.yday ) );
    end
    
    vline(xx2.yday(Itt),'k--')
    
    clear Lot2 Lot1
    Lot1=xx2.Lot(:,Itt);
    Lot2=REsamp.Lot(:,idd,whc);
    
    clear I1 Ip In
    I1=find( isnan(Lot1) & ~isnan(Lot2)  );
    I2=find( ~isnan(Lot1) & isnan(Lot2)  );
    Ip=find(Lot2>Lot1);
    In=find(Lot2<Lot1);
    
    hold on
    timesamp=REsamp.tsamp(:,idd,whc);
    depthsamp=REsamp.zsamp(:,idd,whc);
    plot(timesamp(I1),depthsamp(I1),'mo')
    %     plot(timesamp(I2),depthsamp(I2),'yo')
    %     plot(timesamp(Ip),depthsamp(Ip),'go')
    %     plot(timesamp(In),depthsamp(In),'bo')
end

%%
%if saveplots==1
fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_ExampleFalsePos'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-pdf')
%end

%%

numel(I1)/numel(Lot2)*100
numel(I2)/numel(Lot2)*100
numel(Ip)/numel(Lot2)*100
numel(In)/numel(Lot2)*100

%% contour resampled and corresponding true epsilons to compare, for one realization (start time)

whc=23

figure(3);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 1,3);
cl=[-9 -3]
cmap=jet;
cmap=[ 0.7*[1 1 1];cmap];

idt=isin(xx2.yday,[160 185]);

axes(ax(1))
ezpc(xx2.yday(idt),xx2.z,log10(xx2.eps(:,idt)))
caxis(cl)
%xlim(xl)
colorbar
colormap(cmap)
title(['Tchain ' num2str(whmoor) ', w=' num2str(REsamp.w_samp) 'm/s'],'interpreter','none')
SubplotLetterMW('Tchain')

axes(ax(2))
ezpc(REsamp.timeall(whc,:),REsamp.z,log10(REsamp.eps(:,:,whc)))
caxis(cl)
colorbar
colormap(cmap)
SubplotLetterMW('resamp')

axes(ax(3))
ezpc(REsamp.timeall_true(whc,:),REsamp.z,log10(REsamp.eps_true(:,:,whc)))
caxis(cl)
colorbar
colormap(cmap)
SubplotLetterMW('tchain samp')

linkaxes(ax)

%% Plot mean profile of epsilon w/ and w/o false overturns
% (just for this time period and case number)

e2=REsamp.eps(:,idd,whc);
e2(I1)=nan;

figure(2);clf
agutwocolumn(0.6)
wysiwyg
semilogx(nanmean(xx2.eps(:,Itt),2),REsamp.z)
hold on
semilogx(nanmean(REsamp.eps(:,idd,whc),2),REsamp.z)
semilogx(nanmean(e2,2),REsamp.z)
xlim([1e-9 1e-5])
axis ij
grid on
legend('true','resamp','w/o false')

%%
Ds='stair'
Nm='probability'
Nm='count'
L2=Lot2;L2(I1)=nan;
figure(2);clf
histogram(Lot1,20,'DisplayStyle',Ds,'Normalization',Nm)
hold on
%histogram(Lot1(I2),20,'DisplayStyle',Ds,'Normalization',Nm)
histogram(Lot2(I1),20,'DisplayStyle',Ds,'Normalization',Nm)
histogram(L2(:),20,'DisplayStyle',Ds,'Normalization',Nm)

%histogram(Lot2(In),20,'DisplayStyle',Ds,'Normalization',Nm)
xlabel('L')

%
%
%
%
%% Similar to 1st above, find false overturns etc. But do for all realizations in ensemble
% to get total/average numbers
%29 Jan 2015

clear ; close all

whmoor=3 % mooring #
testnum=1
minOT=50
saveplots=0

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true' T-chain data
%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% empty arrays for indices
FalsePos   =[];
LotFalsePos=[];
FalseNeg   =[];
ExPos      =[];
LotExPos   =[];
ExNeg      =[];
LotExNeg   =[];
LotSamp    =[];
LotReal    =[];

EpsSamp    =[];
EpsReal    =[];

% time-mean profiles
EpsProfReal=nan * ones( length(REsamp.z) , 100) ;
EpsProfSamp=nan * ones( length(REsamp.z) , 100) ;
EpsProfNoFalse=nan * ones( length(REsamp.z) , 100) ;

for whc=1:100
    
    clear Lot2 Lot1
    % true profiles corresponding to times of resampled profiles
    Lot1=REsamp.Lot_true(:,:,whc);
    % resampled profiles
    Lot2=REsamp.Lot(:,:,whc);
    
    clear I1 Ip In I2
    I1=find( isnan(Lot1) & ~isnan(Lot2)  );
    I2=find( ~isnan(Lot1) & isnan(Lot2)  );
    Ip=find(Lot2>Lot1);
    In=find(Lot2<Lot1);
    
    FalsePos=[FalsePos ; I1(:) ];
    FalseNeg=[FalseNeg ; I2(:) ];
    ExPos=[ExPos ; Ip(:) ];
    ExNeg=[ExNeg ; In(:) ];
    
    LotFalsePos=[LotFalsePos ; Lot2(I1)];
    LotExPos=[LotExPos ; Lot2(Ip)];
    LotExNeg=[LotExNeg ; Lot2(In)];
    
    % all resampled Lot
    iov=find(~isnan(Lot2));
    LotSamp=[LotSamp ; Lot2(iov)];
    
    clear iov
    % all true Lot
    iov=find(~isnan(Lot1));
    LotReal=[LotReal ; Lot1(iov)];
    clear iov
    
    EpsProfReal(:,whc)=nanmean(REsamp.eps_true(:,:,whc),2);
    EpsProfSamp(:,whc)=nanmean(REsamp.eps(:,:,whc),2);
    
    clear e2
    e2=REsamp.eps_true(:,:,whc);
    e2(I1)=nan;
    EpsProfNoFalse(:,whc)=nanmean(e2,2);
    
end % case#

%% total # of points
Ntotal=numel(REsamp.Lot(:,:,1))*100

% percent of points that are false positives
Pfalsepos=numel(FalsePos)/Ntotal *100

% percent of false negatives
Pfalseneg=numel(FalseNeg)/Ntotal *100

% percent of points where Lot is over/under estimated
Pexpos=numel(ExPos)/Ntotal*100
Pexneg=numel(ExNeg)/Ntotal*100

% percent of overturns in resampled data that are false (true profile was
% stable)
numel(FalsePos) / numel(LotSamp) * 100
%%
numel(FalsePos)/numel(FalseNeg)

% total # points with overturns
ida=find(~isnan(REsamp.Lot_true));
Nover=numel(ida)
% total # points withOUT overturns
idb=find(isnan(REsamp.Lot_true));
Nstable=numel(idb)

numel(ida)/numel(REsamp.Lot_true)*100
numel(idb)/numel(REsamp.Lot_true)*100

%numel(ida)/numel(idb)
numel(idb)/numel(ida)
%%
Nstable/Nover

numel(ExPos)/Nstable *100
numel(ExNeg)/Nover *100

%%

for whc=1:100
Lot=REsamp.Lot_true(:,:,whc);
N=numel(Lot);
ida=find(isnan(Lot));
numel(ida)/N*100
pause(1)
end
%%
%~ Plot histogram of overturn sizes
Ds='stair'
%Ds='bar'

Nm='count'
Nm='probability'

Nbar=40
figure(1);clf
agutwocolumn(0.6)
wysiwyg
h=histogram(LotReal,Nbar,'DisplayStyle',Ds,'Normalization',Nm)
hold on
histogram(LotSamp,Nbar,'DisplayStyle',Ds,'Normalization',Nm)
histogram(LotFalsePos,Nbar,'DisplayStyle',Ds,'Normalization',Nm)
ylabel('count','fontsize',16)
xlabel('L (m) ','fontsize',16)
%title('False Pos. Lot')
grid on
legend('Real','resamp','falsepos')
title(['Tchain' num2str(whmoor) ', w=' num2str(REsamp.w_samp) 'm/s'])

%
if saveplots==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_LotHistFalsePos'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
end

% also plot true tchain data (from all times)
%idt=isin(xx2.yday,[nanmin(REsamp.timeall(:)) nanmax(REsamp.timeall(:))]);

%~ Plot time-average profiles of epsilon
figure(2);clf
agutwocolumn(0.6)
wysiwyg
semilogx(nanmean(EpsProfReal,2),REsamp.z,'k.-','linewidth',2)
hold on
semilogx(nanmean(EpsProfSamp,2),REsamp.z,'.-','linewidth',2)
semilogx(nanmean(EpsProfNoFalse,2),REsamp.z,'.--','linewidth',2)
%semilogx(nanmean(xx2.eps(:,idt),2),REsamp.z,'linewidth',2)
legend('real','samp','nofalse')
axis ij
grid on
xlim([1e-9 10^(-5.5)])
set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
%xlim([1e-9 1e-6])
ylim([REsamp.z(1) REsamp.z(end)])
title(['Tchain' num2str(whmoor) ', w=' num2str(REsamp.w_samp) 'm/s'])
ylabel('Depth (m)','fontsize',16)
xlabel('\epsilon (Wkg^{-1}) ','fontsize',16)

if saveplots==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_EpsProfFalsePos'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
end


%
%
%
%% Same as above, but for depth-integrated epsilon
% find false overturns etc. But do for all realizations in ensemble
% to get total/average numbers
% 2 Feb 2015

clear ; close all

whmoor=3 % mooring #
testnum=1
minOT=50
saveplots=0

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true' T-chain data
%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

% load resampled data
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% empty arrays for indices
FalsePos   =[];
LotFalsePos=[];
FalseNeg   =[];
ExPos      =[];
LotExPos   =[];
ExNeg      =[];
LotExNeg   =[];
LotSamp    =[];
LotReal    =[];

EpsSamp    =[];
EpsReal    =[];
%

% depth-integrated timeseries of epsilon

EpsTsReal=nan * ones( length(REsamp.timeall(1,:)) , 100) ;
EpsTsSamp=nan * ones( length(REsamp.timeall(1,:)) , 100) ;
EpsTsNoFalse=nan * ones( length(REsamp.timeall(1,:)) , 100) ;

% common time-vector to use
tcommon=REsamp.timeall(1,:);

for whc=1:100
    
    clear Lot2 Lot1
    % true profiles corresponding to times of resampled profiles
    Lot1=REsamp.Lot_true(:,:,whc);
    % resampled profiles
    Lot2=REsamp.Lot(:,:,whc);
    
    clear I1 Ip In I2
    I1=find( isnan(Lot1) & ~isnan(Lot2)  );
    I2=find( ~isnan(Lot1) & isnan(Lot2)  );
    Ip=find(Lot2>Lot1);
    In=find(Lot2<Lot1);
    
    FalsePos=[FalsePos ; I1(:) ];
    FalseNeg=[FalseNeg ; I2(:) ];
    ExPos=[ExPos ; Ip(:) ];
    ExNeg=[ExNeg ; In(:) ];
    
    LotFalsePos=[LotFalsePos ; Lot2(I1)];
    LotExPos=[LotExPos ; Lot2(Ip)];
    LotExNeg=[LotExNeg ; Lot2(In)];
    
    % all resampled Lot
    iov=find(~isnan(Lot2));
    LotSamp=[LotSamp ; Lot2(iov)];
    
    clear iov
    % all true Lot
    iov=find(~isnan(Lot1));
    LotReal=[LotReal ; Lot1(iov)];
    clear iov
    
    %~~~~~~~~~~~~~~~
    % compute depth-integrated epsilon (time-averaged)
    t_avg=2/24;
    dz_samp=nanmean(diff(REsamp.z));
    
    clear Eps1 Esps2
    Eps1=REsamp.eps(:,:,whc);
    Eps2=REsamp.eps_true(:,:,whc);
    
    % depth-integrate epsilon
    clear epsiint_real epsint_samp epsint_nofalse
    epsint_samp=nansum(Eps1)*dz_samp*1026;
    epsint_real=nansum(Eps2)*dz_samp*1026;
    
    Eps2(I1)=nan;
    epsint_nofalse=nansum(Eps2)*dz_samp*1026;
    
    % do a moving average
    clear eps_samp t_samp eps_real t_real
    [eps_samp,t_samp]=SimpleBoxCar(epsint_samp,t_avg,REsamp.timeall(whc,:));
    [eps_real,t_real]=SimpleBoxCar(epsint_real,t_avg,REsamp.timeall_true(whc,:));
    [eps_nofalse,t_nofalse]=SimpleBoxCar(epsint_nofalse,t_avg,REsamp.timeall_true(whc,:));
    
    % interp to same time vector
    EpsTsSamp(:,whc)=interp1(t_samp,eps_samp,tcommon);
    EpsTsReal(:,whc)=interp1(t_real,eps_real,tcommon);    
    EpsTsNoFalse(:,whc)=interp1(t_samp,eps_nofalse,tcommon);
    %~~~~~~~~~~~~~~~
    
end % case#

figure(1);clf
%semilogy(tcommon,EpsTsReal,'.','color',0.7*[1 1 1])
%hold on
semilogy(tcommon,nanmean(EpsTsReal,2),'linewidth',2)
hold on
semilogy(tcommon,nanmean(EpsTsSamp,2),'linewidth',2)
semilogy(tcommon,nanmean(EpsTsNoFalse,2),'--','linewidth',2)
grid on
title(['Tchain' num2str(whmoor) ', w=' num2str(REsamp.w_samp) 'm/s'])
xlabel('Yearday')
ylabel('\int \epsilon dz')
legend('real','samp','nofalse')


%%

figure(3);clf
e1=REsamp.eps(:,:,30);

semilogy(REsamp.timeall(30,:),e1(80,:),'o')

e2=e1(80,:);

%%

figure(1);clf

loglog(EpsTsSamp(:),EpsTsReal(:),'.')
hold on
%loglog(EpsTsNoFalse(:),EpsTsReal(:),'.')
grid on
xvec=linspace(1e-6,1e1,100);
hold on
loglog(xvec,xvec,'k--')
axis equal
%%

Ds='Stair'
Ds='bar'

Nm='Probability'
%Nm='count'

ib=find(log10(EpsTsReal(:))<-4.8);
EpsTsReal(ib)=nan;

ib=find(log10(EpsTsSamp(:))<-4.8);
EpsTsSamp(ib)=nan;

ib=find(log10(EpsTsNoFalse(:))<-4.8);
EpsTsNoFalse(ib)=nan;

figure(1);clf
histogram(log10(EpsTsReal(:)),'DisplayStyle',Ds,'Normalization',Nm,'edgecolor','none')
hold on
histogram(log10(EpsTsSamp(:)),'DisplayStyle',Ds,'Normalization',Nm,'edgecolor','none')
Ds='stair'
histogram(log10(EpsTsNoFalse(:)),'DisplayStyle',Ds,'Normalization',Nm,'edgecolor','k')
ylabel(Nm,'fontsize',16)
xlabel('\rho \int \epsilon dz','fontsize',16)
grid on
xlim([-5 1])

%% scatterplot bias vs true magnitude

ib=find(log10(EpsTsReal(:))<-4.8);
EpsTsReal(ib)=nan;

ib=find(log10(EpsTsSamp(:))<-4.8);
EpsTsSamp(ib)=nan;

bias=(EpsTsSamp(:)-EpsTsReal(:)) ./EpsTsReal(:) *100;
bias2=(EpsTsNoFalse(:)-EpsTsReal(:))./EpsTsReal(:) *100;
E2=(EpsTsNoFalse(:)-EpsTsReal(:));

figure(1);clf
loglog(EpsTsReal(:),(bias),'.')
%hold on
%loglog(EpsTsReal(:),(bias2),'.')
%loglog(EpsTsReal(:),(E2),'.')
hol don
grid on
legend('samp','nofalse')
xlim([1e-5 1e1])
ylim([1e0 1e6])
%axis equal
hline(1e2,'k--')
ylabel('% error')%,' fontsize',16)
%%


Esamp=REsamp.eps;
ia=find(log10(Esamp)<-10);
Esamp(ia)=nan;

Ereal=REsamp.eps_true;
ib=find(log10(Ereal)<-10);
Ereal(ib)=nan;

figure(1);clf
agutwocolumn(0.7)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.08,1,2);

axes(ax(1))

Ds='Stair'
%Ds='bar'
Nm='Probability'
Nm='count'

histogram(log10(Ereal(:)),30,'DisplayStyle',Ds,'Normalization',Nm)%,'edgecolor','none')
hold on
histogram(log10(Esamp(:)),30,'DisplayStyle',Ds,'Normalization',Nm)%,'edgecolor','none')
%Ds='stair'
%histogram(log10(EpsTsNoFalse(:)),'DisplayStyle',Ds,'Normalization',Nm,'edgecolor','k')
ylabel(Nm,'fontsize',16)
xlabel('\rho \int \epsilon dz','fontsize',16)
grid on
xlim([-9.5 -4])

axes(ax(2))

Ds='Stair'
%Ds='bar'
Nm='Probability'
%Nm='count'

histogram(log10(Ereal(:)),30,'DisplayStyle',Ds,'Normalization',Nm)%,'edgecolor','none')
hold on
histogram(log10(Esamp(:)),30,'DisplayStyle',Ds,'Normalization',Nm)%,'edgecolor','none')
%Ds='stair'
%histogram(log10(EpsTsNoFalse(:)),'DisplayStyle',Ds,'Normalization',Nm,'edgecolor','k')
ylabel(Nm,'fontsize',16)
xlabel('\rho \int \epsilon dz','fontsize',16)
grid on
xlim([-9.5 -4])

%%

bias=(Esamp(:)-Ereal(:)) ./ Ereal(:) * 100;

figure(2);clf
agutwocolumn(0.6)
wysiwyg
loglog(Ereal(:),abs(bias),'.')
%loglog(Ereal(:),(bias),'.')
grid on
ylim([1e-1 1e6])

figure(2);clf
histogram(log10((bias(bias>0))),20,'Normalization','probability')
hold on
histogram(log10(abs(bias(bias<0))),20,'Normalization','probability')
grid on
xlim([-1 5])
%%