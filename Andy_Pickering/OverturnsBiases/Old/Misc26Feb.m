%
%
% Compare profiles for a short time period (one resampled profile) to try
% to get a sense of what is going on
%
%
%
%%


clear ; close all

plotiso=1
whmoor=3
minOT=50
testnum=2

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% xl=[168 170]
% xl=[168.2 168.7]
% xl=[169.2 169.7]
xl=[170.2 170.7]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
xl=[181.7 182.2]
xl=[182.2 182.7]
xl=[180 185]
%xl=[183.2 183.7]
%xl=[184.2 184.7]
xl=[210 212]
xl=[210.1 210.5]
%xl=[210.3 210.6]


figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.1, 1,2);
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
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ' - w=' num2str(REsamp.w_samp) ],'interpreter','none')

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:8:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)


% choose one resampling realization out of ensemble
whc=30

% plot an example sampling track over the true data
axes(ax(1))
hold on
plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w')


clear id_thisperiod id_oneP time_thisP id_true_thisP
id_thisperiod=isin(REsamp.timeall(whc,:),xl);

id_oneP=id_thisperiod(6);

plot(REsamp.tsamp(:,id_oneP,whc),REsamp.zsamp(:,id_oneP,whc),'b')
shg
%
time_thisP=[REsamp.tsamp(1,id_oneP,whc) REsamp.tsamp(end,id_oneP,whc)];
freqline(time_thisP(1))
freqline(time_thisP(2))
shg


axes(ax(2))
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14;
caxis(cl)
xlim(time_thisP)
title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' - ' num2str(xl(2)) ' - w=' num2str(REsamp.w_samp) ],'interpreter','none')

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:4:end),'k')
end

hold on
plot(REsamp.tsamp(:,id_oneP,whc),REsamp.zsamp(:,id_oneP,whc),'b')

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)
xlabel('Yearday','fontsize',16)
%~~

id_true_thisP=isin(xx2.yday,time_thisP);

%~~ plot profiles of temp, N2, and eps from this period
figure(2);clf
ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.06, 3,1);

axes(ax(1))
plot(xx2.T(:,id_true_thisP),xx2.z)
hold on
plot(REsamp.t(:,id_oneP),REsamp.z,'linewidth',3,'color',0.0*[1 1 1])
axis ij
grid on
xlim([2 8])
xlabel('T [^o C]')
ylabel('Depth [m]')

axes(ax(2))
plot(xx2.n2(:,id_true_thisP),xx2.z)
hold on
plot(REsamp.n2(:,id_oneP),REsamp.z,'linewidth',3,'color',0.0*[1 1 1])
axis ij
grid on
xlim(1e-5*[-0.5 4])
xlabel('N^2 [s^{-2}]')
ylabel('Depth [m]')

axes(ax(3))
semilogx(xx2.eps(:,id_true_thisP),xx2.z,'.')
hold on
semilogx(nanmean(xx2.eps(:,id_true_thisP),2),xx2.z,'r','linewidth',2)
semilogx(REsamp.eps(:,id_oneP),REsamp.z,'linewidth',2,'color',0.0*[1 1 1])
xlim([1e-11 1e-4])
axis ij
grid on
xlabel('\epsilon [Wkg^{-1}]')
ylabel('Depth [m]')

%%