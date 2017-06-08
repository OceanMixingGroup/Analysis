%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ExamineShortSectionsResamp.m
%
%
% 7 Jan 2015
% 17 Feb 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot entire timeseries first to pick sections

% 19 Dec. 2014

clear ; close all

plotiso=1
minOT=50
whmoor=1 % mooring #

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

% first plot all data to choose sections
figure(1);clf ; agutwocolumn(1) ; wysiwyg

cl=[-9 -4]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

ezpc(xx2.yday,xx2.z,log10(xx2.eps))
cb=colorbar
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14
%
caxis(cl)
title(['T-chain ' num2str(whmoor) ],'fontsize',16)

if plotiso==1
    hold on
    %contour(xx2.yday(id),xx2.z,xx2.T(:,id),[1:7],'k')
    contour(xx2.yday,xx2.z,xx2.T,tm(1:18:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)

%% Plot only a shorter section

clear ; close all

plotiso=0
whmoor=4 % mooring #
minOT=50

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

load( fullfile( 'Data' , ['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps']) )


xl=[168 170]
xl=[168.3 168.6]
%xl=[178.2 178.5]
%xl=[179 179.6]
xl=[179.6 180.2]
xl=[181 183]
%xl=[180 185]
%xl=[183.2 183.7]
%xl=[210 212]


figure(2);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 2,3);
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

axes(ax(1))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14
%
caxis(cl)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ],'interpreter','none')
%title(['T-chain ' num2str(whmoor) ],'fontsize',16)

if plotiso==1
    hold on
    %contour(xx2.yday(id),xx2.z,xx2.T(:,id),[1:7],'k')
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:18:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])

ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)

%~~ load 1 resampled dataset to get sizes etc.
% whmoor=3
testnum=1
whcase=5
clear fname x_resamp
%fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]   ;
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
load(fullfile('Data/',fname))

hold on
plot(x_resamp.tsampall,x_resamp.zall,'w')
%
idr=isin(x_resamp.time,xl);

%idr=idr(2:5)
vline(x_resamp.time(idr),'k--')
shg

figure(3);clf
ax2 = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 2,1);

axes(ax2(1))
plot(x_resamp.T(:,idr),x_resamp.z)
axis ij
grid on
xlabel('T')

axes(ax2(2))
plot(x_resamp.Lot(:,idr),x_resamp.z,'o')
axis ij
grid on
xlabel('L_T')

linkaxes(ax2,'y')

%  load all resampled data to get mean profiles etc.

%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%figure(4);clf
axes(ax(2))
eps_mean=nan*ones(length(REsamp.z),100);
for whc=1:100
    idt=isin(REsamp.timeall(whc,:),xl);
    e2=squeeze(REsamp.eps(:,idt,whc));
    semilogx(nanmean(e2,2),REsamp.z,'.','color',0.75*[1 1 1]);
    hold on
    eps_mean(:,whc)=nanmean(e2,2);
end

axis ij
xlim([1e-9 1e-4])
hr=semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
% plot mean of all resampling
hre=semilogx(nanmean(eps_mean,2),REsamp.z,'linewidth',2)
legend([hr hre],'True','Resamp','location','best')
grid on
%ylabel('Depth (m) ','fontsize',16)
xlabel('\epsilon ','fontsize',16)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ', w_samp=' num2str(REsamp.w_samp) 'm/s' ],'interpreter','none')
%
id2=isin(REsamp.timeall(1,:),xl);
id=isin(xx2.yday,xl);

Lot_All_re=reshape(REsamp.Lot(:,id2,:),numel(REsamp.Lot(:,id2,:)),1);
Lttot_All_re=reshape(REsamp.Lttot(:,id2,:),numel(REsamp.Lot(:,id2,:)),1);
eps_All_re=reshape(REsamp.eps(:,id2,:),numel(REsamp.eps(:,id2,:)),1);

Lot_All_tc=reshape(xx2.Lot(:,id),numel(xx2.eps(:,id)),1);
Lttot_All_tc=reshape(xx2.Lttot(:,id),numel(xx2.eps(:,id)),1);
eps_All_tc=reshape(xx2.eps(:,id),numel(xx2.eps(:,id)),1);

eps_All_re(log10(eps_All_re)<=-11)=nan;
eps_All_tc(log10(eps_All_tc)<=-11)=nan;


Nm='probability';
%Nm='pdf'l
Ds='Stair';

axes(ax(3))
h1=histogram(log10(eps_All_tc(:)),20,'Normalization',Nm,'Displaystyle',Ds)
xlabel('Log_{10} \epsilon','fontsize',16)
hold on
h2=histogram(log10(eps_All_re(:)),20,'Normalization',Nm,'Displaystyle',Ds)
xlabel('Log_{10} \epsilon (Wkg^{-1})','fontsize',16)
grid on
ylabel('Probability')
legend([h1 h2],'true','samp','location','best')

axes(ax(4))
h1=histogram((Lot_All_tc(:)),'Normalization',Nm,'Displaystyle',Ds)
xlabel('L (m)','fontsize',16)
hold on
h2=histogram((Lot_All_re(:)),'Normalization',Nm,'Displaystyle',Ds)
grid on
ylabel('Probability')
legend([h1 h2],'true','samp','location','best')

axes(ax(5))
h1=histogram((Lttot_All_tc(:)),'Normalization',Nm,'Displaystyle',Ds)
xlabel('L_T (m)','fontsize',16)
hold on
h2=histogram((Lttot_All_re(:)),'Normalization',Nm,'Displaystyle',Ds)
grid on
ylabel('Probability')
legend([h1 h2],'true','samp','location','best')

%
% need to get real T-chain data at the times of each resampled profiles

%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_AllCases']))

LotAll_real=[];
LotAll_samp=[];

%figure(5) ; clf

for whc=1:100
    clear idd tvec
    tvec=REsamp.timeall(whc,:);
    idd=isin(tvec,xl);
    
    clear Itt
    Itt=nan*ones(length(idd),1);
    for tt=1:length(idd)
        [val,Itt(tt)]=nanmin( abs( tvec(idd(tt)) - xx2.yday ) );
    end
    
    clear Lot2 Lot1
    Lot2=x_resamp.Lot(:,idd);
    Lot1=xx2.Lot(:,Itt);
    
    LotAll_real=[LotAll_real ; Lot1(:) ];
    LotAll_samp=[LotAll_samp ; Lot2(:) ];
end

axes(ax(6))
histogram(LotAll_samp(:)-LotAll_real(:))
xlim(600*[-1 1])
xlabel('L_{samp}-L_{true}')


%% Different version - 10 Feb

clear ; close all

plotiso=1
whmoor=3 % mooring #
minOT=50

testnum=1

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps']) )
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

xl=[168 170]
xl=[168.2 168.7]
xl=[169.2 169.7]
xl=[170.2 170.7]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
xl=[181.7 182.2]
xl=[182.2 182.7]
%xl=[180 185]
%xl=[183.2 183.7]
%xl=[184.2 184.7]
%xl=[210 212]




figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 1,3);
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

% plot true T-chain data
axes(ax(1))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14
caxis(cl)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ],'interpreter','none')
% plot isopycnals
if plotiso==1
    hold on
    %contour(xx2.yday(id),xx2.z,xx2.T(:,id),[1:7],'k')
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:18:end),'k')
end
cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)


%~~ load resampled dataset
clear fname x_resamp
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

whc=80

%~ Plot the T-chain epsilon sub-sampled with the period of sampling (but
% still instantenous vertical profiles)
axes(ax(2))
ezpc(REsamp.timeall_true(whc,:),REsamp.z,log10(REsamp.eps_true(:,:,whc)))
vline(REsamp.timeall_true(whc,:),'k--')
xlim(xl)
caxis(cl)
colorbar


%~ Plot the resampled epsilon for one case
axes(ax(3))
%idt=isin(REsamp.timeall_true(whc,:),xl);
%ezpc(REsamp.timeall_true(whc,idt),REsamp.z,log10(REsamp.eps_true(:,idt,whc)))
ezpc(REsamp.timeall(whc,:),REsamp.z,log10(REsamp.eps(:,:,whc)))
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')
vline(REsamp.timeall(whc,:),'k--')
xlim(xl)
caxis(cl)
colorbar

% plot the sampling track over the true data
axes(ax(1))
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')

linkaxes(ax)

%% Another Different version - 17 Feb
% contour average ensemble

clear ; close all

plotiso=1
whmoor=3
minOT=50
testnum=1

% xl=[168 170]
% xl=[168.2 168.7]
% xl=[169.2 169.7]
% xl=[170.2 170.7]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
%xl=[181 182.2]
%xl=[182.2 182.7]
xl=[180 185]
%xl=[183.2 183.7]
%xl=[184.2 184.7]
%xl=[210 212]
%xl=[210 210.6]
%xl=[210.25 210.4]
%xl=[195 200]
xl=[205 215]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))



%~~ average epsilon over entire ensemble (all phases)

[M,N]=size(REsamp.timeall);
dt=nanmean(diff(REsamp.timeall(1,:)));

% interp resampled fields from each 'case' to same time
tc=REsamp.timeall(M,1):dt:REsamp.timeall(1,end);
ei=nan*ones(length(REsamp.z),length(tc),M);
Lti=ei;
Loti=ei;

for whc=1:M
    
    for whz=1:length(REsamp.z)
        ei  (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.eps(whz,:,whc)  ,tc);
        Lti (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.Lttot(whz,:,whc),tc);
        Loti(whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.Lot(whz,:,whc)  ,tc);
    end
    
end


eint_samp=squeeze(nansum(ei,1));
eint_real=nansum(xx2.eps);
t_avg=2/24
eint_samp_avg=nanmean(eint_samp,2);
[Esamp,Tsamp]=SimpleBoxCar(eint_samp_avg',t_avg,tc);
[Ereal,Treal]=SimpleBoxCar(eint_real,t_avg,xx2.yday);



figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.01, 1,3);
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

% plot true T-chain data
axes(ax(2))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14;
caxis(cl)
SubplotLetterMW('\epsilon true')

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:18:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xtloff
%xlabel('Yearday','fontsize',16)

whc=10
% plot an example sampling track over the true data
axes(ax(2))
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')
xlim(xl)


%~~
axes(ax(3))
ezpc(tc,REsamp.z,log10(nanmean(ei,3)))
xlim(xl)
caxis(cl)
colorbar
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
SubplotLetterMW('<\epsilon> samp')

% ** Plot timeseries of integrated epsilon also

axes(ax(1))
semilogy(tc,eint_samp,'.','color',0.7*[1 1 1])
hold on
hs=semilogy(Tsamp,Esamp,'linewidth',2)
semilogy(tc,nanmean(eint_samp,2),'o','color',0.5*[1 1 1])
hold on
hr=semilogy(Treal,Ereal,'k')
cb=colorbar;killcolorbar(cb)
ylim([1e-8 10^(-2.5)])
xlim(xl)
grid on
%xlabel('Yearday')
ylabel('\int \epsilon','fontsize',16)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ' - w=' num2str(REsamp.w_samp) ],'interpreter','none')
legend([hs hr],'Sampled','True','location','best')
xtloff

linkaxes(ax,'x')




% plot depth profiles

% * plot indiivudal also in gray to see spread
% * check that mean of resamptrue equals mean of true

idsamp=isin(tc,xl);
figure(2);clf
agutwocolumn(0.7)
wysiwyg
em=nanmean(ei,3);
semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'linewidth',2)
hold on
%semilogx(nanmean(ei(:,idsamp,:),2),REsamp.z)
semilogx(nanmean(em(:,idsamp),2),REsamp.z,'linewidth',2)

%emt=nanmean(REsamp.eps_true(:,idsamp,:),3);
%semilogx(nanmean(emt,2),REsamp.z,'--')

axis ij
xlim([1e-9 1e-5])
grid on
legend('Tchain','Sampled')
xlabel('<\epsilon>','fontsize',16)
ylabel('Depth (m)','fontsize',16)
title(['Tchain ' num2str(whmoor) ' yday ' num2str(xl(1)) ' - ' num2str(xl(2))])

%% Plot histograms for this period

Ncases=size(REsamp.eps,3)

Lot_samp=[];
Lt_samp=[];
eps_samp=[];
Otnsq_samp=[];
for whc=1:Ncases
    clear idd L1  L2 n2
    idd=isin(REsamp.timeall(whc,:),xl)
    L1= REsamp.Lot_each(:,idd,whc);
    Lot_samp=[Lot_samp ; L1(~isnan(L1))];
    
    L2= REsamp.Lt_each(:,idd,whc);
    Lt_samp=[Lt_samp ; L2(~isnan(L2))];
    
    e2= REsamp.eps_each(:,idd,whc);
    eps_samp=[eps_samp ; e2(~isnan(e2))];
    
    n2= REsamp.Otnsq_each(:,idd,whc);
    Otnsq_samp=[Otnsq_samp ; n2(~isnan(n2))];
    
end

%%

% *** 20 fEb 2015 ***

Nm='pdf'
%Nm='count'

figure(8);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 2,2);

axes(ax(1))
histogram(Lot_samp(:),'Normalization',Nm)
hold on
histogram(xx2.Lot_each(:,id),'Normalization',Nm)
histogram(REsamp.Lot_each(:),'Normalization',Nm,'DisplayStyle','stair')
histogram(xx2.Lot_each(:),'Normalization',Nm,'DisplayStyle','stair')
legend('samp','true','sampall','trueall')
xlabel('L_{ot}','fontsize',16)

axes(ax(2))
histogram(Lt_samp(:),'Normalization',Nm)
hold on
histogram(xx2.Lt_each(:,id),'Normalization',Nm)
histogram(REsamp.Lt_each(:),'Normalization',Nm,'DisplayStyle','stair')
histogram(xx2.Lt_each(:),'Normalization',Nm,'DisplayStyle','stair')
xlabel('L_{T}','fontsize',16)

axes(ax(3))
histogram(log10(eps_samp(:)),'Normalization',Nm)
hold on
histogram(log10(xx2.eps_each(:,id)),'Normalization',Nm)
histogram(log10(REsamp.eps_each(:)),'Normalization',Nm,'DisplayStyle','stair')
hold on
histogram(log10(xx2.eps_each(:)),'Normalization',Nm,'DisplayStyle','stair')
xlabel('log_{10}\epsilon','fontsize',16)

axes(ax(4))
histogram(Otnsq_samp(:),'Normalization',Nm)
hold on
histogram(xx2.Otnsq_each(:,id),'Normalization',Nm)
histogram(REsamp.Otnsq_each(:),'Normalization',Nm,'DisplayStyle','stair')
histogram(xx2.Otnsq_each(:),'Normalization',Nm,'DisplayStyle','stair')
xlabel('N^2','fontsize',16)
xlim([0 0.5e-5])

%%



%% 5 Mar. 2015 -  contour average ensemble for different sampling speeds to see effect

clear ; close all

plotiso=1
whmoor=3
minOT=50

% xl=[168 170]
% xl=[168.2 168.7]
% xl=[169.2 169.7]
% xl=[170.2 170.7]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
%xl=[181 182.2]
%xl=[182.2 182.7]
xl=[180 185]
xl=[182 184]
xl=[183.2 183.6]
%xl=[183.2 183.7]
%xl=[184.2 184.7]
%xl=[210 212]
%xl=[210 210.6]
%xl=[210.25 210.4]
xl=[196.2 196.6]
xl=[197.2 197.6]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~~ load resampled dataset
clear fname x_resamp testnum M N dt tc ei idsamp ep
testnum=3
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% average epsilon over entire ensemble (all phases)

[M,N]=size(REsamp.timeall);
dt=nanmean(diff(REsamp.timeall(1,:)));

% interp resampled fields from each 'case' to same time
tc=REsamp.timeall(M,1):dt:REsamp.timeall(1,end);
ei=nan*ones(length(REsamp.z),length(tc),M);
for whc=1:M
    for whz=1:length(REsamp.z)
        ei  (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.eps(whz,:,whc)  ,tc);
    end
end
%~~~

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.01, 1,5);
cl=[-9 -3]

tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

% plot true T-chain data
axes(ax(1))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14;
caxis(cl)
SubplotLetterMW('\epsilon true')
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ],'interpreter','none')
% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:18:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xtloff
%xlabel('Yearday','fontsize',16)

idsamp=isin(tc,xl);
figure(2);clf
semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k')
hold on
clear ep
ep=squeeze(nanmean(ei(:,idsamp,:),2));
semilogx(nanmean(ep,2),REsamp.z)
axis ij
%
whc=10


%~~
axes(ax(2))
ezpc(tc,REsamp.z,log10(nanmean(ei,3)))
xlim(xl)
caxis(cl)
colorbar
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
%SubplotLetterMW('<\epsilon> samp')
SubplotLetterMW(['w=' num2str(REsamp.w_samp) 'm/s'])

% plot an example sampling track over the true data
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')



%~~~ load another resampled dataset
clear fname x_resamp testnum M N dt tc ei idsamp ep
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% average epsilon over entire ensemble (all phases)

[M,N]=size(REsamp.timeall);
dt=nanmean(diff(REsamp.timeall(1,:)));

% interp resampled fields from each 'case' to same time
tc=REsamp.timeall(M,1):dt:REsamp.timeall(1,end);
ei=nan*ones(length(REsamp.z),length(tc),M);
for whc=1:M
    for whz=1:length(REsamp.z)
        ei  (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.eps(whz,:,whc)  ,tc);
    end
end
%~~~


%~~
axes(ax(3))
ezpc(tc,REsamp.z,log10(nanmean(ei,3)))
xlim(xl)
caxis(cl)
colorbar
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
%SubplotLetterMW('<\epsilon> samp')
SubplotLetterMW(['w=' num2str(REsamp.w_samp) 'm/s'])

% plot an example sampling track over the true data
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')


idsamp=isin(tc,xl);
figure(2);
hold on
clear ep
ep=squeeze(nanmean(ei(:,idsamp,:),2));
semilogx(nanmean(ep,2),REsamp.z)
axis ij



%~~~ load another resampled dataset
clear fname x_resamp testnum M N dt tc ei idsamp ep
testnum=2
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% average epsilon over entire ensemble (all phases)

[M,N]=size(REsamp.timeall);
dt=nanmean(diff(REsamp.timeall(1,:)));

% interp resampled fields from each 'case' to same time
tc=REsamp.timeall(M,1):dt:REsamp.timeall(1,end);
ei=nan*ones(length(REsamp.z),length(tc),M);
for whc=1:M
    for whz=1:length(REsamp.z)
        ei  (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.eps(whz,:,whc)  ,tc);
    end
end
%~~~


%~~
axes(ax(4))
ezpc(tc,REsamp.z,log10(nanmean(ei,3)))
xlim(xl)
caxis(cl)
colorbar
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
%SubplotLetterMW('<\epsilon> samp')
SubplotLetterMW(['w=' num2str(REsamp.w_samp) 'm/s'])
%
% plot an example sampling track over the true data
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')

idsamp=isin(tc,xl);
figure(2);
hold on
clear ep
ep=squeeze(nanmean(ei(:,idsamp,:),2));
semilogx(nanmean(ep,2),REsamp.z)
axis ij


%~~~ load another resampled dataset
clear fname x_resamp testnum M N dt tc ei idsamp ep
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% average epsilon over entire ensemble (all phases)

[M,N]=size(REsamp.timeall);
dt=nanmean(diff(REsamp.timeall(1,:)));

% interp resampled fields from each 'case' to same time
tc=REsamp.timeall(M,1):dt:REsamp.timeall(1,end);
ei=nan*ones(length(REsamp.z),length(tc),M);
for whc=1:M
    for whz=1:length(REsamp.z)
        ei  (whz,:,whc)=interp1(REsamp.timeall(whc,:),REsamp.eps(whz,:,whc)  ,tc);
    end
end
%~~~

%~~
axes(ax(5))
ezpc(tc,REsamp.z,log10(nanmean(ei,3)))
xlim(xl)
caxis(cl)
colorbar
xlabel('Yearday','fontsize',16)
ylabel('Depth','fontsize',16)
%SubplotLetterMW('<\epsilon> samp')
SubplotLetterMW(['w=' num2str(REsamp.w_samp) 'm/s'])

% plot an example sampling track over the true data
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')

%
idsamp=isin(tc,xl);
figure(2);
hold on
clear ep
ep=squeeze(nanmean(ei(:,idsamp,:),2));
semilogx(nanmean(ep,2),REsamp.z)
axis ij
legend('true','1','2','3','4')
xlim([1e-9 1e-5])
grid on
xlabel('<\epsilon>','fontsize',16)
ylabel('Depth','fontsize',16)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ],'interpreter','none')
%


linkaxes(ax)

%
%
%%
%%
figure(2);clf
semilogx(nanmean(eps_prof_all_up,2),REsamp.z)
hold on
%semilogx((eps_prof_all_up),REsamp.z)
semilogx(nanmean(eps_prof_all_dn,2),REsamp.z)
%semilogx((eps_prof_all_dn),REsamp.z)
axis ij

idtrue=isin(xx2.yday,xl);
semilogx(nanmean(xx2.eps(:,idtrue),2),xx2.z,'k')
xlim([1e-9 1e-5])
grid on
%ylim([REsamp.z(1) REsamp.z(end)])
ylim([1000 REsamp.z(end)])
ylabel('Depth','fontsize',16)
xlabel('<\epsilon>','fontsize',16)
%%
