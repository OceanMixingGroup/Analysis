%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% GetAllSpecEQ14.m
%
% Script to get all spectra from EQ14 files, so I can try to examine some
% composite spectra from different epsilon ranges etc.
%
% The spectra are saved during the chi-pod calculations.
%
% * also save fits
%
%----------------
% 08/03/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
clear ; close all

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP'

fmax=7

all_k=[];
all_spec=[];

all_k_fit=[];
all_spec_fit=[];

all_eps_cham=[];
all_chi_cham=[];
all_castnum=[];

all_chi_chi=[];
all_eps_chi=[];

dz=50

for castnum=1:3100
    
    disp(['Working on castnum ' num2str(castnum)])
    clear avg avg1
    try
        
        load( fullfile(datdir,['zsm10m_fmax' num2str(fmax) 'Hz_respcorr0_fc_99hz_gamma20'],['EQ14_' sprintf('%04d',castnum) 'avg.mat']) )
        avg1=avg;clear avg
        
        % load chameleon profile also to compare
        
        clear cham idb chameps chamchi
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',castnum) '.mat'])
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',castnum) '.mat'])
        cham=avg;clear avg
        
        % don't use chis where epsilon is NaN
        idb=find(isnan(cham.EPSILON));
        cham.CHI(idb)=nan;
        
        chameps=interp1(cham.P,cham.EPSILON,avg1.P);
        chamchi=interp1(cham.P,cham.CHI,avg1.P);        
        
        all_k=[all_k ; avg1.ks(1:dz:end,:)];
        all_k_fit=[all_k_fit ; avg1.kks(1:dz:end,:)];

        all_spec=[all_spec ; avg1.kspec(1:dz:end,:)];
        all_spec_fit=[all_spec_fit ; avg1.kkspec(1:dz:end,:)];
        
        all_eps_cham=[all_eps_cham ; chameps(1:dz:end)];
        all_chi_cham=[all_chi_cham ; chamchi(1:dz:end)];

        all_chi_chi=[all_chi_chi ; avg1.chi1(1:dz:end)];
        all_eps_chi=[all_eps_chi ; avg1.eps1(1:dz:end)];
                
    end % try
    
end % castnum

% Interpolate all spectra to a common k vector
clear ki speci
ki=linspace(0,300,300);
speci=nan*ones(length(all_eps_cham),length(ki));
speci_fit=speci;%nan*ones(length(all_eps_cham),length(ki));
for whz=1:length(all_eps_cham)
    try
        speci(whz,:)=interp1(all_k(whz,:),all_spec(whz,:),ki);
        speci_fit(whz,:)=interp1(all_k_fit(whz,:),all_spec_fit(whz,:),ki);
    end
end

AllSpec=struct()
AllSpec.all_k=all_k;
AllSpec.all_k_fit=all_k_fit;
AllSpec.all_spec=all_spec;
AllSpec.all_spec_fit=all_spec_fit;
AllSpec.all_eps_cham=all_eps_cham;
AllSpec.all_chi_cham=all_chi_cham;
AllSpec.all_chi_chi=all_chi_chi;
AllSpec.all_eps_chi=all_eps_chi;
AllSpec.ki=ki;
AllSpec.speci=speci;
AllSpec.speci_fit=speci_fit;
save(['AllSpec_' num2str(dz) '.mat'],'AllSpec')


%% Plot all spectra and the mean

figure(1);clf
loglog(AllSpec.ki,AllSpec.speci,'.')
hold on
loglog(AllSpec.ki,nanmedian(AllSpec.speci),'k')

%% Plot for different ranges of epsilon (from chameleon)

ig1=find( log10(AllSpec.all_eps_cham)<-8 );%& log10(AllSpec.all_eps_cham)>-8.2);
ig2=find( log10(AllSpec.all_eps_cham)>-8 );

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(211);
loglog(AllSpec.ki,AllSpec.speci(ig1,:),'.')
hold on
loglog(AllSpec.ki,nanmedian(AllSpec.speci(ig1,:)),'k','linewidth',2)
grid on
ylim([1e-12 1e0])
freqline(2)
freqline(7)
title('log_{10}\epsilon < -8')

ax2=subplot(212);
loglog(AllSpec.ki,AllSpec.speci(ig2,:),'.')
hold on
loglog(AllSpec.ki,nanmedian(AllSpec.speci(ig2,:)),'k','linewidth',2)
grid on
ylim([1e-12 1e0])
freqline(2)
freqline(7)
title('log_{10}\epsilon > -8')
xlabel('wavenumber')

linkaxes([ax1 ax2])

%%
% for whz=1:100:length(all_chi_cham)
%     figure(1);clf
%     loglog(all_k(whz,:),all_spec(whz,:))
%     xlim([1e0 1e3])
%     ylim([1e-12 1e0])
%     pause(0.1)
% end

%%

ig1=find(log10(all_eps_cham)<-8);
ig2=find(log10(all_eps_cham)>-8);

figure(2);clf

subplot(211)
imagesc(log10(speci(ig1,:)))
colorbar
caxis([-11 -0])
xlim([0 60])
title('log_{10}\epsilon < -8')

subplot(212)
imagesc(log10(speci(ig2,:)))
colorbar
caxis([-11 -0])
xlim([0 60])
title('log_{10}\epsilon > -8')

%%


