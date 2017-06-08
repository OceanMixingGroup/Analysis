%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Misc2_9Nov.m
%
% This code loads some data from various resampling scenarios I ran, and
% plots the results as well as the actual full T-chain data, to see if
% there is a bias etc.
%
% Test cases were done in Misc8Nov.m and results saved as mat files that
% are loaded in here. Each test case # has the same sampling shifted in
% time by ~20mins.
%
%
% 10 NOv. - note: need to recompute epsilon for Tchain (gives different
% answer than the eps that is already in file?)
%
% AP 10 Nov. 2014
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%trange=,[x.time(1) x.time(end)];
trange=[164 171]
trange=[177 187] % 'b'
cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

addpath /Users/Andy/Cruises_Research/LADCP_processing/ctd_proc2/


% load in 'real' data
load('all_moorings.mat')
load('all_gridded.mat')
c=[3];
xx=grd{c}
clear grd


% recompute Tchain epsilon (eps field already in that structure is not correct?)
xx.yday=datenum2yday(xx.time);
idtr=isin(xx.yday,trange)

xx2=xx;
xx2.time=xx.time(idtr);
xx2.yday=xx.yday(idtr);
xx2.T=xx.T(:,idtr);
xx2.S=34.604-.045*(xx2.T-2.5);
xx2.eps=NaN*xx2.S;

hb=waitbar(0)
for ind=1:length(xx2.time)
    waitbar(ind/length(xx2.time),hb)
        if mean(xx2.T(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,1,1e-5,0);
        xx2.eps(:,ind)=Epsout;
    end
  
end

delete(hb)
%

%% Plot data for 1 case and true value

%
whcase=2
load(['Tchain3_resamp_case' num2str(whcase) 'b.mat'])
x=x_resamp;
%
figure(31);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.01, 1,4);

axes(ax(1));
ezpc(xx2.yday,xx2.z,xx2.T)
colorbar
caxis([0 8])

axes(ax(2));
ezpc(x.time(1:end-1),x.z,x.T)
colorbar
caxis([0 8])

axes(ax(3));
ezpc(xx2.yday,xx2.z,log10(xx2.eps))
caxis([-9 -2])
colorbar

axes(ax(4));
ezpc(x.time(1:end-1),x.z,log10(x.eps))
caxis([-9 -2])
colorbar

linkaxes(ax)

%% Plot time-averaged depth profiles of epsilon

cols=['m' 'r' 'g' 'b' 'm' 'k' 'r' 'g' 'b' 'm']
figure(2);clf
agutwocolumn(0.75)
wysiwyg

emall=[];

for whcase=1:10
    load(['Tchain3_resamp_case' num2str(whcase) 'b.mat'])
    x=x_resamp;    
    semilogx(nanmean(x.eps,2),x.z,'color',cols(whcase))
    %plot(nanmean(x.eps,2),x.z,'color',cols(whcase))
    emall=[emall nanmean(x.eps,2)];
    hold on
    axis ij
end

% plot mean of all test cases
semilogx(nanmean(emall,2),x.z,'color',0.5*[1 1 1],'linewidth',3)

% plot 'real' profile
semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',3)
xlim([1e-7 1e-5])
ylim([1100 2100])
ylabel('Depth','fontsize',16)
xlabel('<\epsilon> (Wkg^{-1})','fontsize',16)
title(['ydays ' num2str(x.time(1)) ' - ' num2str(x.time(end))],'fontsize',16)

%fname='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/Tchain3_164_171_Profiles'
%save2pdfAP(fname)



%% plot timeseries of depth-average epsilon


cols=['k' 'r' 'g' 'b' 'm' 'k' 'r' 'g' 'b' 'm']
figure(2);clf
agutwocolumn(0.75)
wysiwyg

tavg=0.2
for whcase=1:10
    load(['Tchain3_resamp_case' num2str(whcase) '.mat'])
    x=x_resamp;
    
    [eps_out,time_out]=SimpleBoxCar(nanmean(x.eps),tavg,x.time);
    
    %    semilogx(nanmean(x.eps,2),x.z,'color',cols(whcase))
    semilogy(time_out,eps_out,'color',cols(whcase))
    %    semilogy(x.time(1:end-1),nanmean(x.eps),'color',cols(whcase))
    hold on
    %    axis ij
end

% also plot real eps
% load('all_moorings.mat')
% load('all_gridded.mat')
% c=[3];
% xx=grd{c}
% clear grd

%
t_real=datenum2yday(xx2.time);
z_real=xx2.z;
%idtr=isin(t_real,[x.time(1) x.time(end)]);
[eps_out,time_out]=SimpleBoxCar(nanmean(xx2.eps),tavg,t_real);
semilogy(time_out,eps_out,'k','linewidth',2)
%xlim([x.time(1) x.time(end)])
grid on
xlabel('Yearday','fontsize',16)
ylabel('<\epsilon>','fontsize',16)
title(['ydays ' num2str(x.time(1)) ' - ' num2str(x.time(end))],'fontsize',16)
%%
fname='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/Tchain3_164_171_EpsavgTS'
save2pdfAP(fname)

%
%
%
%% Plot timeseries of depth-integrated epsilon


cols=['k' 'r' 'g' 'b' 'm' 'k' 'r' 'g' 'b' 'm']
figure(2);clf
agutwocolumn(0.75)
wysiwyg

%tavg=0.25
tavg=1
eintall=[];
tintall=[];
for whcase=1:10
    clear x esum eps_out time_out
    load(['Tchain3_resamp_case' num2str(whcase) 'b.mat'])
    x=x_resamp;
    x.time
    esum=    nansum(x.eps)*2;
    [eps_out,time_out]=SimpleBoxCar(esum,tavg,x.time);
    eintall=[eintall ; eps_out];
    tintall=[tintall ; time_out];
    %  semilogy(x.time(1:end-1),esum,'o-','color',cols(whcase))
    semilogy(time_out,eps_out,'-','color',cols(whcase))
    hold on
    %    axis ij
end

%
t_real=datenum2yday(xx2.time);
z_real=xx2.z;
idtr=isin(t_real,[x.time(1) x.time(end)]);
[eps_out,time_out]=SimpleBoxCar(nansum(xx2.eps)*10,tavg,t_real);
semilogy(time_out,eps_out,'k','linewidth',3)

semilogy(tintall(1,:),nanmean(eintall),'color',0.5*[1 1 1],'linewidth',3)

grid on
xlabel('Yearday','fontsize',16)
ylabel('<\int{\epsilon}dz>','fontsize',16)
title(['ydays ' num2str(x.time(1)) ' - ' num2str(x.time(end))],'fontsize',16)
xlim(trange)

%%

fname='/Users/Andy/Cruises_Research/IWISE/Notes/OverturnBiases/Tchain3_164_169_EpsIntTS'
save2pdfAP(fname)

%%

nanmax(eintall)./nanmean(eintall)


%%

figure(3);clf
plot(tintall(1,:),nanmean(eintall),'color',0.5*[1 1 1],'linewidth',3)
hold on
for whcase=1:10
plot(tintall(whcase,:),eintall(whcase,:))
end

%%
figure(55);clf
plot(tintall(1,:),nanstd(eintall)./nanmean(eintall))
%% interpolate real data to same grid and plot error?

figure(7);clf
for whcase=1:10
    
ei= interp1(tintall(whcase,:),eintall(whcase,:),time_out);

perr=100*(ei-eps_out)./eps_out;
%plot(time_out,ei-eps_out,'o-','color',cols(whcase))
plot(time_out,perr,'o-','color',cols(whcase))
hold on
    
end


%%
