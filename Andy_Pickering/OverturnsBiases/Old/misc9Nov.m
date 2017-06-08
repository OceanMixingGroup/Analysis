%%

% work on a better way to resample data (instead of inefficient loop)
%
% take advantage of fact that depth-vector and corresponding indices in
% real data is repeated each profile (every other is flipped)

%%

clear zvec Nz t_samp t_indsmat It
Nprof=20
%~~ vector of sample depths for each profile (will be same for all, every
%other flipped in direction)
zvec=z_range(1) : dz_samp : z_range(end) ; 

%~~ # of points in one profile
Nz=length(zvec);

dt_samp_day=dt_samp/86400;
%~~ make time vector of sampling for all profiles
t_samp=t0 : dt_samp_day : t0 + dt_samp_day*((Nz*Nprof)-1 );

t_indsmat=nan*ones(length(zvec),length(t_samp));
It=nan*ones(size(t_samp));
hb=waitbar(0)
for wht=1:length(t_samp);
    waitbar(wht/length(t_samp),hb)
    [val,It(wht)]=nanmin(abs(t_samp(wht)-t_real));
    
end
delete(hb)
%
t_indsmat=reshape(It,length(zvec),Nprof);


%%

Iz=nan*ones(size(zvec));
for whz=1:length(zvec)
    [val,Iz(whz)]=nanmin(abs(zvec(whz)-z_real))    ;
end
Iz=Iz(:);

z_indsmat=repmat(Iz,1,Nprof);
%%


clear data_resamp
data_resamp=nan*ones(size(z_indsmat));
for whp=1:Nprof
data_resamp(:,whp)=data_real(z_indsmat(:,whp),t_indsmat(:,whp));
end


%% 
Izvec=repmat(Iz,Nprof,1);
data_resamp=nan*ones(size(It))
for ii=1:length(It)
data_resamp(ii)=data_real(Izvec(ii) , It(ii)  );    
end

data_re_mat=reshape(data_resamp,nt_per_prof,Nprof);

%%


idtr=isin(xx.yday,[t_samp(1) t_samp(end)]);

figure
subplot(211)
ezpc(t_real(idtr),z_real,data_real(:,idtr))
caxis([1 8])
xlim(xl)
subplot(212)
%ezpc(Ts,Zs,data_resamp)
ezpc(t_grid(1:end-1),zvec,data_re_mat)

xlim(xl)
caxis([1 8])


%% Try using interp2 instead of just finding closest grid point

tmp=interp2(t_real,z_real,data_real,t_samp,z_samp)

tmp2=reshape(tmp,nt_per_prof,Nprof);

zvec=z_samp(1:nt_per_prof);

figure(3);clf
ezpc(t_grid(1:end-1),zvec,tmp2)

%%
% %%

idtr=isin(xx.yday,[t_samp(1) t_samp(end)]);
% 
 [T,Z]=meshgrid(t_real(idtr),z_real);
% 
 [Ts,Zs]=meshgrid(t_samp,z_samp);
%
 data_resamp=interp2(T,Z,data_real,Ts,Zs);
%% 
% %%
% xl=[215 215.4]
% figure(1);clf
% subplot(211)
% ezpc(T,Z,data_real)
% caxis([1 8])
% xlim(xl)
% subplot(212)
% ezpc(Ts,Zs,data_resamp)
% xlim(xl)
% caxis([1 8])
% %%