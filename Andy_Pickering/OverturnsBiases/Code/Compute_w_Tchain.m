%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_w_Tchain.m
%
% Compute vertical velocity from T-chain data.
%
%
% 13 Jan 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

minOT=50
whmoor=3;

%~~  Load full T-chain data
load(fullfile('Data',['Tchain' num2str(whmoor) ],['Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT)]))

% compute a low-passed temperature field

addpath /Users/Andy/Cruises_Research/CTDstruct/

% add sgth to Tchain
xx2.sgth=sw_pden(xx2.S,xx2.T,sw_pres(xx2.z,20.5)',0);

% compute displacement
xx2=AddIsopycnalDepthToCTD(xx2,xx2.z,1);

dz=xx2.z(2)-xx2.z(1)
% vertical velocity
xx2.w=diffs(xx2.eta)./dz;
% compute a low-passed temperature field
eta_low=MyLowpass(xx2.yday,xx2.eta,4,24);
% low-frequency vertical velocity
w_low=diffs(eta_low)./dz;
%w_low(abs(w_low)>1.5)=nan;
%%

figure(1);clf
%histogram(xx2.w(:))
histogram(w_low(:))

%%
xl=[180 185]

id=isin(xx2.yday,xl);

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.04, 1,5);

axes(ax(1))
ezpc(xx2.yday(id),xx2.z,xx2.T(:,id));
hold on
plot(xx2.yday(id),xx2.d_Iso(1:15:end,id),'k')
caxis([0 8])
xlim(xl)
colorbar
SubplotLetterMW('T')
%

axes(ax(2))
ezpc(xx2.yday(id),xx2.z,xx2.eta(:,id))
hold on
plot(xx2.yday(id),xx2.d_Iso(1:15:end,id),'k')
%contour(xx2.yday(id),xx2.z,xx2.T(:,id),[0:8],'k')
colorbar
colormap(bluered)
caxis(200*[-1 1])
xlim(xl)
SubplotLetterMW('\eta')

axes(ax(3))
ezpc(xx2.yday(id),xx2.z,xx2.w(:,id))
colorbar
caxis(1*[-1 1])
colormap(bluered)
xlim(xl)
SubplotLetterMW('w')

axes(ax(4))
ezpc(xx2.yday(id),xx2.z,eta_low(:,id))
colorbar
colormap(bluered)
caxis(200*[-1 1])
xlim(xl)
SubplotLetterMW('\eta_{low}')

axes(ax(5))
ezpc(xx2.yday(id),xx2.z,w_low(:,id))
colorbar
caxis(1*[-1 1])
colormap(bluered)
xlim(xl)
SubplotLetterMW('w_{low}')

linkaxes(ax)


%%

figure(3);clf
plot(xx2.yday,xx2.w(50,:))
hold on
plot(xx2.yday,w_low(50,:))
ylim(2*[-1 1])
xlim(xl)
grid on
%%
