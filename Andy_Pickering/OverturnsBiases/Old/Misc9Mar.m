

%%

clear ; close all


load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3/Tchain3_RecomputedEps_MinOT_50.mat')


% make sythetic field

t=0 : 2/3600/24 : 1 ;
%z=0:10:2000;
z=xx2.z;

Tsyn=12 % period of oscillations (hr)
om=2*pi*24/Tsyn;

tmean=nanmean(xx2.T,2);
z=z(~isnan(tmean));
tmean=tmean(~isnan(tmean));

[T,Z]=meshgrid(t,z);

temp=nan*ones(size(T));
temp=repmat(tmean,1,length(t))+ 0.4.*sin(om.*T);

% make a simulated sampling path
t0=0;
z_range=[nanmin(z) nanmax(z)];
dz_samp=10;
time_range=1;
w_samp=0.07;
plotit=0
[tmat,zmat]=MakeProfPath(t0,z_range,dz_samp,time_range,w_samp,plotit);

Nprof=size(zmat,2)

data_resamp=nan*ones(length(z),Nprof);
t_grid=nan*ones(1,Nprof);
%for whp=1:Nprof
whp=2
clear tvec idt datai
tvec=tmat(:,whp);
idt=isin(t,[tvec(1) tvec(end)]);
datai= interp2(t(idt),z,temp(:,idt),tvec,zmat(:,1));

datai2=interp2(t(idt),z,temp(:,idt),tvec,zmat(:,whp));


%%

whp=2
[X,Y] = meshgrid(t(idt),z);
[Xq,Yq] = meshgrid(tvec,zmat(:,whp));
V=temp(:,idt);
Vq = interp2(X,Y,V,Xq,Yq);

[Xq2,Yq2] = meshgrid(tvec,zmat(:,1));
Vq2 = interp2(X,Y,V,Xq2,Yq2);
%%

figure(3);clf
%ezpc(X,Y,V)
ezpc(Xq,Yq,Vq)
%%
datai-diag(Vq)
%%
diag(Vq)
diag(Vq2)
figure(1);clf
plot(diag(Vq),z,flipud(diag(Vq2)),z)
%%
whp=2
figure(1);clf
ezpc(t(idt),z,temp(:,idt))
hold on
contour(T,Z,temp,tmean(1:10:end),'k')
shg
plot(tvec,zmat(:,1),'w')
plot(tvec,zmat(:,whp),'m')
%%
datai-flipud(datai2)

figure(1);clf
plot(datai,zmat(:,1))
hold on
plot(datai2,zmat(:,whp))
axis ij
%%
%datai= interp2(t(idt),z,temp(:,idt),tvec,zmat(:,whp));
%data_resamp(:,whp)=datai(:);
%t_grid(whp)=nanmean(tvec);

%
figure(1);clf
subplot(121)
ezpc(T,Z,temp)
hold on
contour(T,Z,temp,tmean(1:10:end),'k')
cb=colorbar
cb.Label.String='temp'
xlabel('Time (days)')
ylabel('Depth (m) ')

plot(tmat(:,whp),zmat(:,whp),'w')
shg

%
%figure(2);clf
subplot(122)
plot(temp(:,round(nanmean(idt))),z)
hold on
plot(datai,z)
plot(datai,zmat(:,whp))
axis ij
shg
legend('true','samp','samp2')
%%

%%