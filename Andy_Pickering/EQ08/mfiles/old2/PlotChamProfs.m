%%

figure(1);clf

subplot(141)
plot(cal2.fspd,cal2.P)
axis ij
grid on
xlabel('fspd')

subplot(142)
plot(cal2.S,cal2.P)
axis ij
grid on
xlabel('S')

subplot(143)
plot(cal2.T,cal2.P)
axis ij
grid on
xlabel('T')

subplot(144)
plot(cal2.TP,cal2.P)
axis ij
grid on
xlabel('TP')

%%