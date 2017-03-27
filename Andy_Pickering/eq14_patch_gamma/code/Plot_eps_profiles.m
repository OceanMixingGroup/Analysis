%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_eps_profiles.m
%
% Plot profiles of epsilon from chameleon and chi-pod method. Compare
% time-average profiles for groups of casts?
%
% 3/27/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%whN2dTdz = 'line'
whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1

minR2 = 0.0
eq14_patches_paths
dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches')

%% Plot a single profile

cnum=2011

% patch N^2,dTdz w/ constant gamma
load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
ch=avg;clear avg
%avg_patchN2dTdz_constGam=avg;clear avg

% % patch N^2,dTdz w/ patch gamma
% load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
% avg_patchN2dTdzGam=avg;clear avg

% % regular chi-pod method on binned data
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
chb=avg;clear avg

load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])


figure(1);clf
plot(log10(avg.EPSILON),avg.P,'k','linewidth',2)
hold on
plot(log10(chb.eps1),chb.P,'.-','markersize',10,'color',0.5*[1 1 1])
plot(log10(ch.eps1),ch.P,'rp','linewidth',2,'markersize',12)
axis ij
grid on
ylim([0 200])

%% Try plotting average of a group of profiles?

% first find profiles near each other
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

%
figure(1);clf
subplot(211)
plot(cham.castnumber,cham.lon)
subplot(212)
plot(cham.castnumber,cham.lat)

%%

cnums = [1000:1500]
cnums = [2000:3100]

e1=[];
P1=[];
e2=[];
P2=[];
for i=1:length(cnums)
    try
        cnum=cnums(i);
        clear avg
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        e1 = [e1(:) ; avg.eps1(:)];
        P1 = [P1(:) ; avg.P(:)   ];
        
        %load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        e2 = [e2(:) ; avg.eps1(:)];
        P2 = [P2(:) ; avg.P(:)   ];

    end
end

%

iCham=find(cham.castnumber>cnums(1) & cham.castnumber<nanmax(cnums));

%

[dataout1 zout1 Nobs] = binprofile(e1,P1, 0, 10, 200,1);
[dataout2 zout2 Nobs] = binprofile(e2,P2, 0, 10, 200,1);
[cham_bin zout_cham Nobs] = binprofile(cham.EPSILON(:,iCham),cham.P(:,iCham), 0, 10, 200,1);
%%
figure(1);clf
h1 = plot(log10(dataout1),zout1,'bo-','linewidth',2)
hold on
h2 = plot(log10(dataout2),zout2,'ms-','linewidth',2)
h3 = plot(log10(cham_bin),zout_cham,'rp-','linewidth',2)
axis ij
grid on
%xlim([-11 -4])
ylim([0 200])
xlabel('log_{10}[\epsilon]')
legend([h1 h2 h3],'patch','bin','cham','location','best')

eq14_patches_paths
print( fullfile(fig_dir,['eps_prof_example']),'-dpng')
%%