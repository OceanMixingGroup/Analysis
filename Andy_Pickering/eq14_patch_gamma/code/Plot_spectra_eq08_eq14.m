%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_spectra_eq08_eq14.m
%
%
%
%------------------
% 5/17/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%project = 'eq14'
project = 'eq08'

cnum = 1513

eval([project '_patches_paths'])

Params.fmax=32
load( fullfile(path_chipod_bin,['zsm10m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma20_nfft_128'],[project '_' sprintf('%04d',cnum) '_avg.mat']) )
avg1=avg; clear avg

Params.fmax=15
load( fullfile(path_chipod_bin,['zsm10m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma20_nfft_128'],[project '_' sprintf('%04d',cnum) '_avg.mat']) )
avg2=avg; clear avg

Params.fmax=7
load( fullfile(path_chipod_bin,['zsm10m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma20_nfft_128'],[project '_' sprintf('%04d',cnum) '_avg.mat']) )
avg3=avg; clear avg

%%

figure(1);clf
plot(avg1.fstop,avg1.P,'o')
axis ij
hold on
plot(avg2.fstop,avg2.P,'d')
plot(avg3.fstop,avg3.P,'p')

%%

ChkMkDir(fullfile(fig_dir,'spectra'))

%cols = ['b' 'r' 'k']
for iz = 1:10:length(avg1.P)
    
    figure(1);clf
    agutwocolumn(1)
    wysiwyg
    
    subplot(211)
    loglog(avg1.fspec,avg1.tpspec(iz,:),'k')
    grid on
    freqline(32,'b--')
    freqline(15,'r--')
    freqline(7,'m--')
    ylabel('\Phi')
    xlabel('f')
    xlim([1e0 1e2])
    %ylim([1e-7 1e-2])
    
    %figure(2);clf
    subplot(212)
    loglog(avg1.ks(iz,:),avg1.kspec(iz,:),'kd-')
    hold on
    h1 = loglog(avg1.kks(iz,:),avg1.kkspec(iz,:),'b');
    h2 = loglog(avg2.kks(iz,:),avg2.kkspec(iz,:),'r');
    h3 = loglog(avg3.kks(iz,:),avg3.kkspec(iz,:),'m');
    
    ylim([1e-10 1e0])
    xlim([10^(-0.25) 10^(2.5)])
    grid on
    legend([h1 h2 h3],'32hz','15hz','7hz')
    
    % ** multiply by fspd to get k....
    
    fspd = avg1.fspd(iz);
    
    freqline(32/fspd,'b--')
    freqline(15/fspd,'r--')
    freqline(7/fspd,'m--')
    ylabel('\Phi')
    xlabel('k')
    
    title([project ' cast ' num2str(cnum) ' , P=' num2str(avg1.P(iz))])
    
    pause(0.1)
    
    print( fullfile( fig_dir,'spectra',[project_short '_cnum_' num2str(cnum) '_iz_' num2str(iz)]), '-dpng')
    
end



%%

figure(3);clf
plot(avg1.fstop,avg1.P)
hold on
plot(avg2.fstop,avg2.P)
axis ij

%% Try averaging profiles w/ similar chi or epsilon?

id = find(log10(avg1.chi1)>-6 & log10(avg1.chi1)<-4) ;

% interpolate to common k vectors


%%
figure(1);clf
loglog(avg1.fspec,avg1.tpspec(id,:))
hold on
loglog(avg1.fspec,nanmean(avg1.tpspec(id,:)),'k','linewidth',2)


%%