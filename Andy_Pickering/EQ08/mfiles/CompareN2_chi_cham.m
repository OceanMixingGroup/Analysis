%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareN2_chi_cham.m
%
% Compare N2 and dTdz computed in cali_eq08.m to what I compute for
% CTD-chipod processing
%
%
%--------------
% 03/18/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all
Params.z_smooth=1
Ncasts=900
% folder where processed cham files are saved
path_save='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/casts'

for icast=16:Ncasts
    
    try
        
        %~~ load calibrated cast file
        clear fn cal cal2
        fn=['eq08' sprintf('%04d',icast) '.mat']
        load(fullfile(path_save,fn))
        %~~ compute 2nd N2 and dTdz as I would for CTD-chipod
        % average temp and sal in 1m bins like we normally do for CTD data
        clear zmin dz zmax tbin zbin sbin
        zmin=nanmin(cal2.P);
        dz=1;
        zmax=nanmax(cal2.P);
        minobs=2;
        [tbin zbin Nobs] = binprofile(cal2.T ,cal2.P, zmin, dz, zmax,minobs);
        [sbin zbin Nobs] = binprofile(cal2.S,cal2.P, zmin, dz, zmax,minobs);
        %clear zmin dz zmax minobs
        
        % compute dT/dz and N^2 for chi calculations
        clear ctd z_smooth
        ctd=struct();
        ctd.t1=tbin;
        ctd.s1=sbin;
        ctd.p=zbin;
        ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        ctd=Compute_N2_dTdz_forChi(ctd,Params.z_smooth);
        
        [N23 zbin Nobs] = binprofile(cal2.N2 ,cal2.P, zmin, dz, zmax,minobs);
        
        %~~ compare the results
        %%
        figure(1);clf
        semilogx(cal2.N2,cal2.P)
        hold on
        semilogx(ctd.N2,ctd.p)
        semilogx(N23,ctd.p)
        grid on
        axis ij
        
        figure(2);clf
        loglog(ctd.N2,N23,'.')
        hold on
        xvec=linspace(1e-6,1e-3,100);
        loglog(xvec,xvec,'k--')
        grid on
        
        %%
        pause(1)
        
    end % try
    
end % icast

%%

figure(2);clf
loglog(ctd.N2,N23,'.')
hold on
xvec=linspace(1e-6,1e-3,100);
loglog(xvec,xvec,'k--')
grid on

%%
