%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotEpsContribVsL.m
%
%
%
% 10 Feb 2015 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot cumulative epsilon as function of increasing Lot

clear ; close all

saveplot=1
whmoor=3
minOT=50
    
cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

figure(1);clf
agutwocolumn(0.6)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.06, 0.06, 0.1, 0.06, 2,1);

% load 'True' Tchain data
load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

clear idt Lvec B I epsvec eps_sort CS Nvec
%find data in same time period as resampled data
idt=isin(xx2.yday,[180 185]);

% sort overturns
Lvec=xx2.Lot(:,idt);Lvec=Lvec(:);
%Lvec=xx2.Lttot(:,idt);Lvec=Lvec(:);
[B,I] = sort(Lvec,1,'ascend');

% sort eps using same indices
epsvec=xx2.eps(:,idt);epsvec=epsvec(:);
eps_sort=epsvec(I);

% compute cumulative sum
CS=nancumsum(epsvec(I));

% plot vs overturn size
axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)
xlabel('L')
ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon')
grid on
hold on

% do same for Thorpe scale
clear  Lvec B I epsvec eps_sort CS Nvec

% sort overturns
%Lvec=xx2.Lot(:,idt);Lvec=Lvec(:);
Lvec=xx2.Lttot(:,idt);Lvec=Lvec(:);
[B,I] = sort(Lvec,1,'ascend');

% sort eps using same indices
epsvec=xx2.eps(:,idt);epsvec=epsvec(:);
eps_sort=epsvec(I);

% compute cumulative sum
CS=nancumsum(epsvec(I));

% plot vs Thorpe scale
axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)
grid on
hold on

% load Resampled data
clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

% %~~~
% Lvec=REsamp.Lot_true(:);
% %Lvec=REsamp.Lttot_true(:);
% [B,I] = sort(Lvec);
% epsvec=REsamp.eps_true(:);
% CS=nancumsum(epsvec(I));
% 
% axes(ax(1))
% plot(B,CS./nansum(epsvec),'-','linewidth',2)
% hold on
% 
% %Lvec=REsamp.Lot(:);
% Lvec=REsamp.Lttot_true(:);
% [B,I] = sort(Lvec);
% epsvec=REsamp.eps_true(:);
% CS=nancumsum(epsvec(I));
% 
% axes(ax(2))
% plot(B,CS./nansum(epsvec),'-','linewidth',2)
% hold on
% %
%~~~
clear Lvec B I epsvec eps_sort CS
Lvec=REsamp.Lot(:);
%Lvec=REsamp.Lttot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));

axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% do for Thorpe scale
clear Lvec B I epsvec CS
%Lvec=REsamp.Lot(:);
Lvec=REsamp.Lttot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));

axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% do again for a different speeds
clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
Lvec=REsamp.Lot(:);
%Lvec=REsamp.Lttot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));

axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% plot vs Thorpe scale
clear Lvec B I epsvec CS
%Lvec=REsamp.Lot(:);
Lvec=REsamp.Lttot(:);
[B,I] = sort(Lvec);
epsvec=REsamp.eps(:);
CS=nancumsum(epsvec(I));

axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

axes(ax(1))
xlabel('L','fontsize',16)
ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon','fontsize',16)
grid on
title(['T-chain ' num2str(whmoor)])
legend('true','0.25m/s','0.75m/s','location','best')

axes(ax(2))
xlabel('L_T','fontsize',16)
grid on
%
if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_EpsContribVsLot'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
end
%
%
%
%

%% Do same as above, but only for 1 depth range

clear ; close all

saveplot=0
whmoor=4
minOT=50

zrange=[1500 1600]
zrange=[1300 1800]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

figure(1);clf
agutwocolumn(0.6)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.06, 0.06, 0.1, 0.06, 2,1);


clear idt Lvec B I epsvec eps_sort CS Nvec
% find data in same time period as resampled data
%idt=isin(xx2.yday,[180 185]);

% sort overturns
%Lvec=xx2.Lot(:,idt);Lvec=Lvec(:);
%Lvec=xx2.Lttot(:,idt);Lvec=Lvec(:);
% [B,I] = sort(Lvec,1,'ascend');
% 
% % sort eps using same indices
% epsvec=xx2.eps(:,idt);epsvec=epsvec(:);
% eps_sort=epsvec(I);
% 
% % compute cumulative sum
% CS=nancumsum(epsvec(I));
% 
% % plot vs overturn size
% axes(ax(1))
% plot(B,CS./nansum(epsvec),'-','linewidth',2)
% xlabel('L')
% ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon')
% grid on
% hold on
% 
% % do same for Thorpe scale
% clear  Lvec B I epsvec eps_sort CS Nvec
% 
% % sort overturns
% %Lvec=xx2.Lot(:,idt);Lvec=Lvec(:);
% Lvec=xx2.Lttot(:,idt);Lvec=Lvec(:);
% [B,I] = sort(Lvec,1,'ascend');
% 
% % sort eps using same indices
% epsvec=xx2.eps(:,idt);epsvec=epsvec(:);
% eps_sort=epsvec(I);
% 
% % compute cumulative sum
% CS=nancumsum(epsvec(I));
% 
% % plot vs Thorpe scale
% axes(ax(2))
% plot(B,CS./nansum(epsvec),'-','linewidth',2)
% grid on
% hold on

% load Resampled data
clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=1
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

idz=isin(REsamp.z,zrange);
%
%~~~
Lvec=REsamp.Lot_true(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps_true(idz,:);epsvec=epsvec(:);
%Lvec=REsamp.Lttot_true(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));

axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)
hold on

%Lvec=REsamp.Lot_true(idz,:);Lvec=Lvec(:);
Lvec=REsamp.Lttot_true(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps_true(idz,:);epsvec=epsvec(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));

axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)
hold on
%
%~~~
clear Lvec B I epsvec eps_sort CS
Lvec=REsamp.Lot(idz,:);Lvec=Lvec(:);
%Lvec=REsamp.Lttot(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps(idz,:);epsvec=epsvec(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));


axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% do for Thorpe scale
clear Lvec B I epsvec CS
clear Lvec B I epsvec eps_sort CS
%Lvec=REsamp.Lot(idz,:);Lvec=Lvec(:);
Lvec=REsamp.Lttot(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps(idz,:);epsvec=epsvec(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));


axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% do again for a different speeds
clear testnum REsamp Lvec B I epsvec eps_sort CS
testnum=4
load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

clear Lvec B I epsvec eps_sort CS
Lvec=REsamp.Lot(idz,:);Lvec=Lvec(:);
%Lvec=REsamp.Lttot(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps(idz,:);epsvec=epsvec(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));


axes(ax(1))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

% plot vs Thorpe scale
clear Lvec B I epsvec eps_sort CS
%Lvec=REsamp.Lot(idz,:);Lvec=Lvec(:);
Lvec=REsamp.Lttot(idz,:);Lvec=Lvec(:);
epsvec=REsamp.eps(idz,:);epsvec=epsvec(:);
[B,I] = sort(Lvec);
CS=nancumsum(epsvec(I));


axes(ax(2))
plot(B,CS./nansum(epsvec),'-','linewidth',2)

axes(ax(1))
xlabel('L')
ylabel('\Sigma_{0}^{i}\epsilon / \Sigma \epsilon')
grid on
title(['T-chain ' num2str(whmoor)])
legend('true','0.25m/s','0.75m/s','location','best')

axes(ax(2))
xlabel('L_T')
grid on
%
if saveplot==1
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_EpsContribVsLot'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
end

%%