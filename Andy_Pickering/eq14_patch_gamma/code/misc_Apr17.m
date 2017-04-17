%%

%% plot average of groups of profiles, binned

clear ; close all

Params.gamma = 0.2;
Params.fmax=7
Params.z_smooth=10

dz=10; % bin size

eq14_patches_paths

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

cnum_range = [ 1500 2000]


clear cnums
cnums = [cnum_range(1) : cnum_range(2) ];

eps_chi= [];
P_chi  = [];
Tz_chi = [] ;
N2_chi = [];

clear Ngood
Ngood=0;
cnum_good=[];
icham=[];
for i=1:length(cnums)
    try
        cnum=cnums(i);
        
        % binned chipod profile
        clear avg
        load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],['EQ14_' sprintf('%04d',cnum) '_avg.mat']))
        eps_chi = [eps_chi(:) ; avg.eps1(:)];
        P_chi   = [P_chi(:) ; avg.P(:)   ];
        Tz_chi  = [Tz_chi(:) ; avg.dTdz(:) ];
        N2_chi  = [N2_chi(:) ; avg.N2(:) ];
        Ngood = Ngood+1;
        
        cnum_good = [cnum_good(:) ; cnum];
        icham = [ icham  find(cham.castnumber==cnum)];
        
    catch
    end
    
end
%

clear ib
ib = find(log10(cham.EPSILON)<-8.5);
cham.EPSILON(ib) = nan;

[chi_bin zout_chi Nobs] = binprofile(eps_chi,P_chi, 0, dz, 200,1);
[cham_bin zout_cham Nobs] = binprofile(cham.EPSILON(:,icham),cham.P(:,icham), 0, dz, 200,1);



%%

ib = find(log10(eps_chi)>-4);

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(121);
plot(log10(eps_chi),P_chi,'.','color',0.7*[1 1 1])
hold on
plot(log10(eps_chi(ib)),P_chi(ib),'rp')
plot(log10(chi_bin),zout_chi,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])
xlabel('log_{10}\epsilon')
title('chi-pod')

ax2 = subplot(122);
plot(log10(cham.EPSILON(:,icham)),cham.P(:,icham),'.','color',0.7*[1 1 1])
hold on
plot(log10(cham_bin),zout_cham,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])
xlabel('log_{10}\epsilon')
title('cham')

linkaxes([ax1 ax2])


figure(3);clf
subplot(121)
h_chi=plot(log10(chi_bin),zout_chi,'k','linewidth',2)
hold on
h_cham = plot(log10(cham_bin),zout_cham,'m','linewidth',2)
axis ij
grid on
legend([h_chi h_cham],'chi','cham')
title('before screening')
xlim([-10 -3])

%%

% try screening out some spikes in chipod data that give huge epsilons
clear ib
ib=find( medfilt1(Tz_chi,5) ./ Tz_chi  >2 ) ;
eps_chi(ib)=nan;

clear ib
ib = find(log10(N2_chi)>-2.5);
eps_chi(ib)=nan;


[chi_bin zout_chi Nobs] = binprofile(eps_chi,P_chi, 0, dz, 200,1);
[cham_bin zout_cham Nobs] = binprofile(cham.EPSILON(:,icham),cham.P(:,icham), 0, dz, 200,1);


figure(2);clf

ax1 = subplot(121);
plot(log10(eps_chi),P_chi,'.','color',0.7*[1 1 1])
hold on
plot(log10(eps_chi(ib)),P_chi(ib),'rp')
plot(log10(chi_bin),zout_chi,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])

ax2 = subplot(122);
plot(log10(cham.EPSILON(:,icham)),cham.P(:,icham),'.','color',0.7*[1 1 1])
hold on
plot(log10(cham_bin),zout_cham,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])

linkaxes([ax1 ax2])


figure(3);

subplot(122)
h_chi=plot(log10(chi_bin),zout_chi,'k','linewidth',2)
hold on
h_cham = plot(log10(cham_bin),zout_cham,'m','linewidth',2)
axis ij
grid on
legend([h_chi h_cham],'chi','cham')
title('after screening')
xlim([-9 -3])

%%

ib = find(log10(eps_chi)>-4);

figure(4);clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(121);
plot(log10(eps_chi),P_chi,'.','color',0.7*[1 1 1])
hold on
plot(log10(eps_chi(ib)),P_chi(ib),'rp')
plot(log10(chi_bin),zout_chi,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])
xlabel('log_{10}\epsilon')
title('chi-pod')

ax2 = subplot(122);
plot(log10(cham.EPSILON(:,icham)),cham.P(:,icham),'.','color',0.7*[1 1 1])
hold on
plot(log10(cham_bin),zout_cham,'k','linewidth',2)
axis ij
grid on
ylim([0 250])
xlim([-13 1])
xlabel('log_{10}\epsilon')
title('cham')

linkaxes([ax1 ax2])

%%