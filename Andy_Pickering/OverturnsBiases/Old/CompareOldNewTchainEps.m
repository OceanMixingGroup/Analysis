%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareOldNewTchainsEps.m
%
% This was a script comparing overturns already computed in data JN sent me
% to re-computed overturns. They were different, I use re-computed
% overturns for all analysis now.
%
%---------------
% 11/10/2014 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

disp('loading data to resample')
%~~ use JN's T-chains during IWISE
load('all_moorings.mat')
load('all_gridded.mat')
c=[3];
xx=grd{c}
clear grd

%%
xx2=xx;
xx2.S=34.604-.045*(xx.T-2.5);
xx2.eps=NaN*xx2.S;

hb=waitbar(0)
for ind=1:length(xx.time)
    waitbar(ind/length(xx.time),hb)
    if mean(xx.T(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,1,1e-5,0);
        xx2.eps(:,ind)=Epsout;
    end
    
end

%%

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,2);

xx.yday=datenum2yday(xx.time);
axes(ax(1));
ezpc(xx.yday,xx.z,log10(xx.eps))
colorbar
caxis([-8 -4])
xlim([190 192])

axes(ax(2));
ezpc(xx.yday,xx2.z,log10(xx2.eps))
colorbar
caxis([-8 -4])
cmap=jet;
cmap=[0.7*[1 1 1];cmap;]
xlim([190 192])
colormap(cmap)

linkaxes(ax)

%%