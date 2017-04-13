%%

% 2762 bad!

clear ; close all

for cnum = 450%2800%:2900
    
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chipod_method_bin/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) '_avg.mat'])

ib = find(log10(avg.eps1)>-4);

figure(1);clf
agutwocolumn(1)

wysiwyg

clear ib2
ib2=find( medfilt1(avg.dTdz,6) ./ avg.dTdz  >2);

ax1=subplot(221);
plot(log10(avg.eps1),avg.P,'-')
hold on
plot(log10(avg.eps1(ib)),avg.P(ib),'o')
plot(log10(avg.eps1(ib2)),avg.P(ib2),'d')
axis ij
grid on
xlabel('\epsilon')
ylim([0 250])
xlim([-12 -2])

ax2=subplot(222);
plot(log10(avg.dTdz),avg.P,'-')
hold on
plot(log10(avg.dTdz(ib)),avg.P(ib),'o')
%plot(log10(medfilt1(avg.dTdz,5)),avg.P,'p')
plot(log10(avg.dTdz(ib2)),avg.P(ib2),'d')
axis ij
grid on
xlabel('T_z')
xlim([-6 0])

ax3=subplot(223);
plot(log10(avg.N2),avg.P,'-')
hold on
plot(log10(avg.N2(ib)),avg.P(ib),'o')
axis ij
grid on
xlabel('N^2')
xlim([-7 -2])

ax4=subplot(224);
plot(log10(avg.chi1),avg.P,'-')
hold on
plot(log10(avg.chi1(ib)),avg.P(ib),'o')
axis ij
grid on
xlabel('\chi')

linkaxes([ax1 ax2 ax3 ax4],'y')

pause

end
%%

figure(2);clf
boxplot(avg.dTdz ./medfilt1(avg.dTdz,5))
