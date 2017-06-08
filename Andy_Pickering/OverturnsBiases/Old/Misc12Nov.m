
%%

clear ; close all

load Tchain3_day194_196resamp

%x_resamp.eps
hist(log10(x_resamp.eps(:)),50)
%hist(log10(x_resamp.Lot(:)))
%hist(log10(x_resamp.Lttot(:)))
%xlim([-8 -4])

%log10(nanmean(x_resamp.eps(:)))
%%

trange=[x_resamp.time(1) x_resamp.time(end)]
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
xx2.Lot=NaN*xx2.S;
xx2.Lttot=NaN*xx2.S;

hb=waitbar(0)
for ind=1:length(xx2.time)
    waitbar(ind/length(xx2.time),hb)
        if mean(xx2.T(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx2.z',xx2.T(:,ind),xx2.S(:,ind),35.8,0,1,1e-5,0);
        xx2.eps(:,ind)=Epsout;
        xx2.Lot(:,ind)=Lot;
        xx2.Lttot(:,ind)=Lttot;
    end
  
end

delete(hb)

%%
%clear 
figure(1);clf
em=nanmean(xx2.eps,2)
semilogx( em,xx.z)
hold on
semilogx(nanmean(x_resamp.eps,2),x_resamp.z,'m')
axis ij

%%

figure(1);clf

subplot(211)
[N,X]=hist(xx2.Lot(:),50)
bar(X,N)
xlim([0 650])

subplot(212)
hist(x_resamp.Lot(:),X)
%bar(X,
xlim([0 650])
%%