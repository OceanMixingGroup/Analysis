%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_Cham_Profs_AP.m
%
% Combine chi-pod profiles from Chameleon casts into one structure.
%
% Processed files made w/ Calc_Chi_eq08_AP.m
%
%----------------
% 03/15/16 - A.Pickering
% 04/28/16 - AP - Add Params.fc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.z_smooth = 1;   % distance to smooth N^2 and dTdz over
Params.nfft     = 128; % # points to use for spectra in each window
Params.fmax     = 7;   % Max frequency to integrate up to
Params.TPthresh = 1e-6;% Minimum threshold for TP variance
Params.resp_corr= 0;   % correct TP spectra for freq response of thermistor
Params.fc       = 99 ;
Params.gamma    = 0.2  % mixing efficiency

% use default fc=99 for no correction (to make file paths same)
if Params.resp_corr==0
    Params.fc=99;
end

% directory where processed casts are saved
datdirsave=fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/',...
        ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft)] )

%    ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_gamma' num2str(Params.gamma*100) '_nfft' num2str(Params.nfft)]);

Ncasts=100;

% parameters for binning
clear zmin dz zmax tbin zbin sbin
zmin=0;
dz=5;
zmax=300;
minobs=1;
zbins=zmin:dz:zmax ;
%
EmpMat=nan*ones(length(zbins),Ncasts);
TPvar=EmpMat;
chi1=EmpMat;
eps=EmpMat;
KT=EmpMat;
dTdz=EmpMat;
N2=EmpMat;
cnums=nan*(ones(1,Ncasts));

hb=waitbar(0);

for icast=1:Ncasts
    waitbar(icast/Ncasts,hb)
    
    try
        clear avg
        load( fullfile(datdirsave,['eq08_' sprintf('%04d',icast) '_avg.mat' ]) )
        
        [chi1(:,icast)  zbin Nobs] = binprofile(avg.chi1  , avg.P, zmin, dz, zmax,minobs);
        [KT(:,icast)    zbin Nobs] = binprofile(avg.KT1   , avg.P, zmin, dz, zmax,minobs);
        [eps(:,icast)   zbin Nobs] = binprofile(avg.eps1  , avg.P, zmin, dz, zmax,minobs);
        [TPvar(:,icast) zbin Nobs] = binprofile(avg.TP1var, avg.P, zmin, dz, zmax,minobs);
        [dTdz(:,icast)  zbin Nobs] = binprofile(avg.dTdz  , avg.P, zmin, dz, zmax,minobs);
        [N2(:,icast)    zbin Nobs] = binprofile(avg.N2    , avg.P, zmin, dz, zmax,minobs);
        cnums(icast)=icast;
        
        
    end % try
    
end % icast
delete(hb)


%
C     = struct();
C.p   = zbins;
C.eps = eps;
C.chi = chi1;
C.KT  = KT;
C.TPvar = TPvar;
C.dTdz= dTdz;
C.N2  = N2;
C.castnumber = cnums;
C.Params = Params;
C.BinParams.dz   = dz;
C.BinParams.zmin = zmin;
C.BinParams.zmax = zmax;
C.MakeInfo = ['Made ' datestr(now) ' w/ Combine_Cham_Profs_AP.m']
C.Info = ['averaged in ' num2str(dz) ' m bins']

save(fullfile('/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/'...
    ,['chi_all_zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100)]),'C');

%%

figure(1);clf
ezpc(C.castnumber,C.p,log10(C.eps))
caxis([-11 -4])
colorbar

%%