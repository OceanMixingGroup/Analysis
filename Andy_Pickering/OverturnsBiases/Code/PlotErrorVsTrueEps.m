%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotErrorVsTrueEps.m
%
% This code plots the error in the resampled depth-integrated epsilon (averaged over a
% specified time period) versus the magnitude of the true epsilon
%
% Modified from a section of PlotResampResults12Nov.m
%
% 22 Dec. 2014
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot depth-integrated  epsilon for 1 resample set (1 speed)

clear ; close all

saveplot=0

whmoor=3    % which mooring to use
t_avg=2/24  % time to average over (hrs)

testnum=2   % which test number (each on is different speed)

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/

whcase=1
fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)]
clear epsint_samp x_resamp eps_samp t_samp
load( fullfile( 'Data' , fname) )

% Load 'true' T-chain data to resample
load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps']) )
% find indices for same time period as resampled data
idtr=isin(xx2.yday,[x_resamp.time(1) x_resamp.time(end)]);
% Make empty arrays
eall_samp=[];
tall_samp=[];

% get data for each case (different phases)
for whcase=1:100
    
    clear fname epsint_samp x_resamp eps_samp t_samp
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase)];
    load( fullfile( 'Data' , fname) )
    
    dz_samp=nanmean(diff(x_resamp.z));
    % depth-integrate epsilon
    epsint_samp=nansum(x_resamp.eps)*dz_samp*1026;
    % do a moving average
    [eps_samp,t_samp]=SimpleBoxCar(epsint_samp,t_avg,x_resamp.time);
    
    if whcase==1
        t_common=t_samp;
    end
    
    % interp to same time vector
    eps_samp=interp1(t_samp,eps_samp,t_common);
    
    % save results with others
    eall_samp=[eall_samp; eps_samp];
    tall_samp=[tall_samp ; t_common];
    
end

% depth-integrate and average true values
epsint=nansum(xx2.eps)*10*1026;
[eps_true,t]=SimpleBoxCar(epsint,t_avg,xx2.yday);

% t,eps  ,  t_samp, ,nanmean(eall_samp)
% interp to same times

% average over ensemble
ens_mean_eps=nanmean(eall_samp);
% interp to same time vector as true data
esamp2=interp1(t_common,ens_mean_eps,t);
% compute bias
bias=esamp2-eps_true;
%bias=(esamp2-eps)./eps*100;

ip=find(bias>0);
in=find(bias<0);

% plot timeseries of epsilon and bias

figure(1) ; clf
ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.04, 1,2);

axes(ax(1))
plot(t,bias,'o')
ylim(0.5*[-1 1])
grid on
gridxy
ylabel('Bias')

axes(ax(2))
semilogy(t,eps_true,'o-',t,esamp2,'s-')
legend('true','resamp')
grid on
xlabel('Yearday')
ylabel('\int \epsilon dz')

linkaxes(ax,'x')

% make a scatter plot of the error vs magnitude of true epsilon

figure(2) ; clf
loglog(eps_true(ip),abs(bias(ip)),'b.','linewidth',2,'markersize',14)
hold on
loglog(eps_true(in),abs(bias(in)),'bo','linewidth',2,'markersize',5)
xlim([1e-4 1e1])
ylim([1e-4 1e1])
grid on
xlabel('\int \epsilon_{true} dz','fontsize',18)
ylabel('\int \epsilon_{resamp}dz- \int \epsilon_{true}dz','fontsize',18)
hold on
xvec=linspace(1e-4,nanmax(eps_true),150);
loglog(xvec,xvec,'k--')
title(['Tchain ' num2str(whmoor) ' - \int \epsilon dz - ' num2str(t_avg*24) ' hr avg' ])
legend('bias>0','bias<0','location','west')
SubplotLetter(['w_{samp}=' num2str(x_resamp.w_samp) 'm/s'])

if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    fname=fullfile(figdir,['Tchain' num2str(whmoor) '_Testnum' num2str(testnum) '_IntEps_ErrorvsTrue' ])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf','-r200')
end

%
%% Plot depth-integrated epsilon for 4 different sampling speeds (2X2 plot)

% Updated 5 Feb 2015

clear ; close all

saveplot=1
whmoor=4    % which mooring to use
t_avg=2/24  % time to average over (days)
minOT=50;

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/

figure(2) ; clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.05, 0.06, 0.1, 0.04, 2,2);
whax=1;

for testnum=[3 1 2 4]   % which test number (each on is different speed)
    
    % load 1 resample case just to get time limits
    clear fname whcase idtr eall_samp tall_samp
    whcase=1;
    fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)];
    clear epsint_samp x_resamp eps_samp t_samp
    load( fullfile( 'Data' , fname) )
    
    % Load 'true' T-chain data to resample
    load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
    
    % find indices for same time period as the resampled data
    idtr=isin(xx2.yday,[x_resamp.time(1) x_resamp.time(end)]);
    
    % Make empty arrays
    eall_samp=[];
    tall_samp=[];
    
    % get data for each case (different phases)
    for whcase=1:100
        
        clear fname epsint_samp x_resamp eps_samp t_samp
        fname=['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)];
        load( fullfile( 'Data' , fname) )
        
        dz_samp=nanmean(diff(x_resamp.z));
        % depth-integrate epsilon
        epsint_samp=nansum(x_resamp.eps)*dz_samp*1026;
        % do a moving average
        [eps_samp,t_samp]=SimpleBoxCar(epsint_samp,t_avg,x_resamp.time);
        
        if whcase==1
            t_common=t_samp;
        end
        
        % interp to same time vector
        eps_samp=interp1(t_samp,eps_samp,t_common);
        
        % save results with others
        eall_samp=[eall_samp; eps_samp];
        tall_samp=[tall_samp ; t_common];
        
    end
    
    % depth-integrate and average true values
    clear epsint eps t
    epsint_true=nansum(xx2.eps)*10*1026;
    [eps_true,t_true]=SimpleBoxCar(epsint_true,t_avg,xx2.yday);
    
    % average over ensemble
    clear ens_mean_eps esamp2 bias ip in
    ens_mean_eps=nanmean(eall_samp);
    
    % interp to same time vector as true data
    esamp2=interp1(t_common,ens_mean_eps,t_true);
    
    % compute bias
    bias=esamp2-eps_true;
    %bias=(esamp2-eps)./eps*100;
    
    clear ig
    ig=~isnan(bias);
    bias=bias(ig);
    eps_true=eps_true(ig);
    
    ip=find(bias>0);
    in=find(bias<0);
    
    % make a scatter plot of the error vs magnitude of true epsilon
    axes(ax(whax))
    loglog(eps_true(ip),abs(bias(ip)),'.','linewidth',2,'markersize',14)
    hold on
    loglog(eps_true(in),abs(bias(in)),'o','linewidth',2,'markersize',5)
    xlim([1e-5 1e1])
    ylim([1e-5 1e1])
    grid on
    
    if whax>2
        xlabel('\rho \int \epsilon_{o} dz','fontsize',18)
    end
    if ~iseven(whax)
        ylabel('\rho \int \epsilon_{s}dz- \rho \int \epsilon_{o}dz','fontsize',18)
    end
    hold on
    
    % plot a 1:1 line for reference
    %   xvec=linspace(1e-4,nanmax(eps),150);
    %   loglog(xvec,xvec,'k--')
    
    clear  X Y P Yfit
    X=eps_true;Y=abs(bias);
    P=polyfit(log10(X),log10(Y),1)
    Yfit=10^(P(2))*X.^(P(1));
    loglog(X,Yfit,'k.')
    
    clear  X Y Pp Yfit
    X=eps_true;Y=abs(bias);
    Pp=polyfit(log10(X(ip)),log10(Y(ip)),1)
    Yfit=10^(Pp(2))*X.^(Pp(1)) ;
    loglog(X,Yfit,'.','color', [0    0.4470    0.7410])
    
    clear  X Y Pn Yfit
    X=eps_true;Y=abs(bias);
    Pn=polyfit(log10(X(in)),log10(Y(in)),1)
    Yfit=10^(Pn(2))*X.^(Pn(1));
    loglog(X,Yfit,'.','color',[0.8500    0.3250    0.0980])
    
    if whax==1
        title(['Tchain ' num2str(whmoor) ' - ' num2str(t_avg*24) ' hr avg' ])
        legend('bias>0','bias<0','location','southeast')
    end
    
    ht=SubplotLetter(['w_{samp}=' num2str(x_resamp.w_samp) 'm/s']);
    ht.FontSize=14;

    whax=whax+1;
    
end

linkaxes(ax)

if saveplot==1
    figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
    fname=fullfile(figdir,['Tchain' num2str(whmoor) '_4Speeds_IntEps_ErrorvsTrue' ])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
end


%%
figure(33);clf
loglog(eps_true,abs(bias),'.','linewidth',2,'markersize',14)
hold on
%    loglog(eps(in),abs(bias(in)),'o','linewidth',2,'markersize',5)
xlim([1e-4 1e1])
ylim([1e-4 1e1])
grid on

ig=~isnan(bias);
X=eps_true(ig);
Y=abs(bias(ig));
P=polyfit(log10(X),log10(Y),1);
Y2=polyval(P,log10(X));

Yfit=10^(P(2))*X.^(P(1));
loglog(X,Yfit,'r.')
%%

figure(3);clf
loglog(eps_true,esamp2,'.','linewidth',2,'markersize',13)
grid on
hold on
axis tight
xvec=linspace(1e-4,nanmax(eps_true),150);
loglog(xvec,xvec,'k--')
xlabel('\epsilon true')
ylabel('\epsilon resamp')

%%

figure(44);clf
histogram(log10(bias(ip)),20)%,'DisplayStyle','stair')
hold on
histogram(log10(abs(bias(in))),20)%,'DisplayStyle','stair')

%%
figure
histogram((bias),40)%,'DisplayStyle','stair')
hold on
%histogram((abs(bias(in))))%,'DisplayStyle','stair')


%%