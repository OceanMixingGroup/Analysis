%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotBiasVsEps.m
%
% Plot bias/error in resampled overturns versus true overturn size or
% epsilon to see what relationship is.
%
% Shows that the bias (sampled eps - true eps), depth-integrated, is
% proportional to the magnitude of the true epsilon (depth-int)
%
%-----------------
% 11/26/14 - AP - apickering@apl.washington.edu
% 07/27/15 - updating to new data paths and forrmats
% 11/23/15 - AP - updated
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot depth-integrated  epsilon for 1 resample set (1 speed)

clear ; close all

% time period to average over
t_avg=2/24

minOT=50;
% which mooring to use
whmoor=3
% which test # (sampling speed)
testnum=3


eall_re=[];
tall_re=[];

% load resampled data
clear epsint_re x_resamp eps_re t_re fname
datadir=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/',['Tchain' num2str(whmoor)],['Test' num2str(testnum)])
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases.mat'];
load(fullfile(datadir,fname))


% load ' true ' data
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain' num2str(whmoor) '/Tchain' num2str(whmoor) '_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
idtr=isin(xx2.yday,[nanmin(REsamp.tsamp(:)) nanmax(REsamp.tsamp(:))]);
zvec=REsamp.z;

% integrate real epsilon
epsint=nansum(xx2.eps)*10*1026;
[e_int_true,t]=SimpleBoxCar(epsint(:,idtr),t_avg,xx2.yday(idtr));

%  
clear eall_re tall_re
% resampled epsilon (integrated)
eall_re=[];
% time vector of resampled epsilons
%tall_re=[];
Nshift=size(REsamp.eps,3)
dz_samp=nanmean(diff(REsamp.z));

% compute depth-integrated epsilon for each ensemble
for whcase=1:Nshift
    clear epsint_re x_resamp eps_re t_re
    epsint_re=nansum(REsamp.eps(:,:,whcase))*dz_samp*1026;
    [eps_re,t_re]=SimpleBoxCar(epsint_re,t_avg,REsamp.tgrid(whcase,:) ) ;    
    % interp to same time as 'true' data
    eall_re=[eall_re; interp1(t_re,eps_re,t)];
end

%%

e_int_samp=nanmean(eall_re);
% plot the true timeseries
figure(1);clf
semilogy(t,e_int_true,'+-','linewidth',2)
hold on
%semilogy(t_re,nanmean(eall_re),'s-','linewidth',2)
semilogy(t,e_int_samp,'s-','linewidth',2)
grid on
title(['Tchain ' num2str(whmoor) ' Depth-integrated \epsilon, ' num2str(t_avg) ' day avg. '])
xlabel('Yearday','fontsize',16)
ylabel('\int \epsilon \rho dz','fontsize',16)
legend('true','samp','location','best')

%%

figure(2);clf
loglog(e_int_true,e_int_samp,'.')
grid on 
xlabel('\epsilon true')
xlabel('\epsilon sampled')
xvec=linspace(1e-4,nanmax(e_int_true),100);
hold on
loglog(xvec,xvec,'k--')

%% Now make a scatter plot of 'bias' vs true value for depth-integrated epsilon

figure(3);clf

%bias=abs(ei2-eps);
bias=abs(e_int_samp-e_int_true);
%plot(eps,bias,'o','linewidth',2)
loglog(e_int_true,bias,'o','linewidth',2)
%semilogx(eps,100*bias./eps,'o','linewidth',2)
grid on
xlabel('\int\epsilon dz','fontsize',16)
ylabel('bias','fontsize',16)
xvec=linspace(1e-4,nanmax(e_int_true),100);
hold on
loglog(xvec,xvec,'k--')
title([num2str(t_avg*24) ' hr avg'],'fontsize',16)
xlim([1e-6 1e2])
ylim([1e-6 1e2])
%%