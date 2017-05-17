%
%
%
%
%%

clear ; close all


load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chipod_method_bin/zsm10m_fmax32Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_1550_avg.mat')
avg1=avg; clear avg

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chipod_method_bin/zsm10m_fmax15Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_1550_avg.mat')
avg2=avg; clear avg

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chipod_method_bin/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_1550_avg.mat')
avg3=avg; clear avg

%%

figure(1);clf
plot(avg1.fstop,avg1.P,'o')
axis ij
hold on
plot(avg2.fstop,avg2.P,'d')
plot(avg3.fstop,avg3.P,'p')

%%

%cols = ['b' 'r' 'k']
for iz = 1:length(avg1.P)

figure(2);clf
loglog(avg1.ks(iz,:),avg1.kspec(iz,:),'kd-')
hold on
h1 = loglog(avg1.kks(iz,:),avg1.kkspec(iz,:),'b');
h2 = loglog(avg2.kks(iz,:),avg2.kkspec(iz,:),'r');
h3 = loglog(avg3.kks(iz,:),avg3.kkspec(iz,:),'m');

ylim([1e-10 1e0])
xlim([1e-1 1e3])
grid on
legend([h1 h2 h3],'32hz','15hz','7hz')

% ** multiply by fspd to get k....

fspd = avg1.fspd(iz);

freqline(32/fspd,'b--')
freqline(15/fspd,'r--')
freqline(7/fspd,'m--')
ylabel('\Phi')
xlabel('k')

pause

title([])

end

%%

figure(3);clf
plot(avg1.fstop,avg1.P)
hold on
plot(avg2.fstop,avg2.P)
axis ij

%%