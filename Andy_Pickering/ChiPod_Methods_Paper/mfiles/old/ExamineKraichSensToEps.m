%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
% Plot Kraichnan spectra for same epsilon and different epsilons. I want to
% check that the shape at resolved wavenumbers is not sensitve to epsilon;
% which I think is why we can get chi right but have epsilon be so off.
%
%
%---------------------
% 06/03/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%%
% Plot Kraichnan spectra for a range of epsilons, along with average FP07
% spectra, to show that the thermistor rolls off at frequencies less than Kb

figure(1);clf
agutwocolumn(0.7)
wysiwyg

xlim([1e-1 1e4])
ylim([1e-8 1e-2])
xlabel('Frequency [Hz]','fontsize',18)
ylabel('\Phi_{T_z} K^2s^{-2}/s^{-1}','fontsize',18)

nu=1e-6;
b_freq=(10.^(-2:.1:3.5))';
fspd=1
tdif=1.5e-7
chi=4e-10
q=7
hh=[]
for eps=[1e-10 1e-9 1e-8]%  1e-8 1e-6  1e-4]%[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4]
clear kb b_spec h
kb = (((eps./(nu.^3)).^.25 )/2/pi).*sqrt(nu./tdif);
b_spec= kraichnan(nu,b_freq/fspd,kb,tdif,chi,q)/fspd;

h=loglog(b_freq,b_spec)
h2=freqline(kb);set(h2(2),'Color',h.Color);
%loglog(b_freq/fspd,b_spec)
hh=[hh h]
hold on

end

%hl=legend([hh ],'\epsilon=1e-10','\epsilon=1e-8','\epsilon=1e-6','\epsilon=1e-4','location','northeast')
%hl.FontSize=15

freqline(1)
freqline(7)


xlabel('Frequency [Hz]','fontsize',18)
ylabel('\Phi_{T_t} [K^2s^{-2} Hz^{-1}]','fontsize',18)
grid on
set(gca,'FontSize',15)

ylim([1e-12 1e0])
grid on
xlabel('wavenumber [cpm]','fontsize',16)
ylabel('\Phi_{dT/dz}[K^2m^{-2}cpm^{-1}]','fontsize',16)
xlim([0.1 1.5e3])

%%