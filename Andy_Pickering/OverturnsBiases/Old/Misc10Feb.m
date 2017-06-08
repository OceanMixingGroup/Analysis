%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Misc10Feb.m
%
% 10 Feb 2015 - AP 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=4
minOT=50

clear REsamp
testnum=3
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

%Nm='count'
Nm='pdf'

Ds='stair'
%Ds='bar'
%
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.06, 0.06, 0.1, 0.06, 2,3);

axes(ax(1))
histogram(REsamp.Lot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(REsamp.Lot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L','fontsize',15)
title('pdf')

Nm='count'
axes(ax(2))
histogram(REsamp.Lot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(REsamp.Lot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L','fontsize',15)
title('count')
legend('true','samp')

vline(nanmean(REsamp.Lot_true(:)),'b--')
vline(nanmean(REsamp.Lot(:)),'r--')


Nm='pdf'
axes(ax(3))
histogram(REsamp.Lttot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(REsamp.Lttot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L_T','fontsize',15)

Nm='count'
axes(ax(4))
histogram(REsamp.Lttot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(REsamp.Lttot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L_T','fontsize',15)

vline(nanmean(REsamp.Lttot_true(:)),'b--')
vline(nanmean(REsamp.Lttot(:)),'r--')

% vline(nanmedian(REsamp.Lttot_true(:)),'b-')
% vline(nanmedian(REsamp.Lttot(:)),'r-')


igt=find(log10(REsamp.eps_true)>-10);
igs=find(log10(REsamp.eps)>-10);

Nm='pdf'
axes(ax(5))
histogram(log10(REsamp.eps_true(igt)),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(log10(REsamp.eps(igs)),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('log_{10} \epsilon','fontsize',15)

Nm='count'
axes(ax(6))
histogram(log10(REsamp.eps_true(igt)),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(log10(REsamp.eps(igs)),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('log_{10} \epsilon','fontsize',15)

vline(log10(nanmean(REsamp.eps_true(igt))),'b--')
vline(log10(nanmean(REsamp.eps(igs))),'r--')

% vline(log10(nanmedian(REsamp.eps_true(igt))),'b-')
% vline(log10(nanmedian(REsamp.eps(igs))),'r-')


%% compare histograms at diff speeds

~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

whmoor=3
minOT=50

clear REsamp
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
RE1=REsamp;clear REsamp

clear REsamp
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
RE2=REsamp;clear REsamp

%Nm='count'
Nm='pdf'

Ds='stair'
%Ds='bar'
%
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.06, 0.06, 0.1, 0.06, 2,3);

axes(ax(1))
%histogram(REsamp.Lot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
histogram(RE1.Lot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(RE2.Lot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L','fontsize',15)
title('pdf')

Nm='count'
axes(ax(2))
histogram(RE1.Lot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(RE2.Lot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L','fontsize',15)
title('count')
legend('true','samp')

vline(nanmean(RE1.Lot_true(:)),'b--')
vline(nanmean(RE2.Lot(:)),'r--')


Nm='pdf'
axes(ax(3))
histogram(RE1.Lttot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(RE2.Lttot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L_T','fontsize',15)

Nm='count'
axes(ax(4))
histogram(RE1.Lttot_true(:),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(RE2.Lttot(:),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('L_T','fontsize',15)

vline(nanmean(RE1.Lttot_true(:)),'b--')
vline(nanmean(RE2.Lttot(:)),'r--')

% vline(nanmedian(REsamp.Lttot_true(:)),'b-')
% vline(nanmedian(REsamp.Lttot(:)),'r-')


igt=find(log10(RE1.eps_true)>-10);
igs=find(log10(RE2.eps)>-10);

Nm='pdf'
axes(ax(5))
histogram(log10(RE1.eps_true(igt)),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(log10(RE2.eps(igs)),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('log_{10} \epsilon','fontsize',15)

Nm='count'
axes(ax(6))
histogram(log10(RE1.eps_true(igt)),20,'Normalization',Nm,'DisplayStyle',Ds)
hold on
histogram(log10(RE2.eps(igs)),20,'Normalization',Nm,'DisplayStyle',Ds)
xlabel('log_{10} \epsilon','fontsize',15)

vline(log10(nanmean(RE1.eps_true(igt))),'b--')
vline(log10(nanmean(RE2.eps(igs))),'r--')

% vline(log10(nanmedian(REsamp.eps_true(igt))),'b-')
% vline(log10(nanmedian(REsamp.eps(igs))),'r-')


%%

for whc=1:100

clear Lot tsamp zsamp ig
Lot=REsamp.Lot(:,:,whc);
tsamp=REsamp.tsamp(:,:,whc);
zsamp=REsamp.zsamp(:,:,whc);

ig=find(Lot>700);

figure(1);clf
ezpc(REsamp.timeall(whc,:),REsamp.z,Lot)
hold on
plot(tsamp(ig),zsamp(ig),'ko')
colorbar
caxis([0 700])

pause(1)
end

%% See if average profile changes much if we exclude the small
% number of overturns larger than ~700m

prof=nan*ones(length(REsamp.z),100);
prof2=prof;
proftrue=prof;

for whc=1:100
    prof(:,whc)=nanmean(REsamp.eps(:,:,whc),2);
    clear ig eps
    ig=find(REsamp.Lot(:,:,whc)>680);
    eps=REsamp.eps(:,:,whc);
    eps(ig)=1e-11;
    prof2(:,whc)=nanmean(eps,2);
    proftrue(:,whc)=nanmean(REsamp.eps_true(:,:,whc),2);
end

figure(1);clf
semilogx(prof,REsamp.z,'color',0.7*[1 1 1])
hold on
semilogx(prof2,REsamp.z,'color',0.3*[1 1 1])
axis ij

figure(2);clf
semilogx(nanmean(proftrue,2),REsamp.z)
hold on
semilogx(nanmean(prof,2),REsamp.z)
semilogx(nanmean(prof2,2),REsamp.z)
xlim([1e-9 1e-6])
grid on
axis ij

%%

[val,iz]=nanmin(abs(REsamp.z-1600))

figure(1);clf
histogram(REsamp.Lot(iz,:,:))
hold on
histogram(REsamp.Lot_true(iz,:,:))

figure(2);clf
histogram(REsamp.Lot(iz,:,:))
hold on
histogram(REsamp.Lot_true(iz,:,:))
%%

%%