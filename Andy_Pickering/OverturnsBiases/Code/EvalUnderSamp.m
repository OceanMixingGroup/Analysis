%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% EvalUnderSamp.m
%
% Evaluate effect of undersampling on overturns
% by sampling real data (instanteous profiles) at
% time intervals of a simulated profiling instrumnet.
%
%------------------
% 25 Nov. 2014 - AP - apickering@apl.washington.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot profiles of time-average epsilon for under-sampling at different speeds

clear ; close all

whmoor=3 ; xl=[1e-9 10^(-5.5)]
whmoor=4 ; xl=[1e-9 10^(-5.9)]
minOT=50
saveplot=1

% load Tchain data
load(fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))

% dt of Tchain data
dt=xx2.yday(2)-xx2.yday(1);

idtr=isin(xx2.yday,[165 185]);

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 2,2);

allws=[0.15 0.25 0.5 0.75];

for whax=1:4 % do for 4 different profiling speeds
    
    clear w_samp Tsamp Tsamp_day nt Ntsamp
    
    % specify vertical sampling velocity
    w_samp=allws(whax);
    
    % time for 1 profile
    Tsamp=range(xx2.z)/w_samp; % sec
    Tsamp_day=Tsamp/86400; % days
    nt=round(Tsamp_day/dt); % ~ # real profiles in 1 resamp profile
    
    % total number of profiles
    Ntsamp=floor(length(xx2.yday(idtr))/nt)-1;
    %    Ntsamp=floor(length(xx2.yday)/nt)-1;
    
    axes(ax(whax))
    semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',2)
    hold on
    axis ij
    grid on
    if ~iseven(whax)
        ylabel('Depth','fontsize',16)
    else
        ytloff
    end
    if whax>2
        xlabel('<\epsilon>','fontsize',16)
    end
    
    % Repeat a bunch of times, shifting start time by 2 mins to get range of possible results
    eall=nan*ones(length(xx2.z),nt);
    for qt=1:nt
        clear idre
        idre=(1+qt): nt : Ntsamp*nt;
        if length(idre)==Ntsamp ; idre=idre(1:end-1) ; end % make all same size
        semilogx(nanmean(xx2.eps(:,idre),2),xx2.z,'.','color',0.7*[1 1 1],'linewidth',1)
        eall(:,qt)=nanmean(xx2.eps(:,idre),2);
        length(idre);
    end
    
    semilogx(nanmean(xx2.eps(:,idtr),2),xx2.z,'k','linewidth',2)
    %    hold on
    % plot mean of re/under-sampled profiles
    semilogx(nanmean(eall,2),xx2.z,'y--','linewidth',1)
    
    % plot +/- standard devation also
    semilogx(nanmean(eall,2)-1.0*std(eall,0,2),xx2.z,'--','color',0.3*[1 1 1])
    semilogx(nanmean(eall,2)+1.0*std(eall,0,2),xx2.z,'--','color',0.3*[1 1 1])
    
    xlim(xl)
    ylim([xx2.z(1) xx2.z(end)])
    ht=SubplotLetter(['w=' num2str(w_samp) ' m/s'],0.6,0.1)
    set(ht,'fontsize',14)
    set(gca,'xtick',[1e-8 1e-7 1e-6 1e-5])
    
end

shg
%
if saveplot==1
    %
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_UnderSamp_4speeds'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    export_fig(fname,'-pdf')
    %
end

%% Same as above, but do for TIMESERIES of depth-integrated epsilon

clear ; close all

whmoor=4
saveplot=0

% load Tchain data
load(fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','Data',['Tchain' num2str(whmoor) '_RecomputedEps']))

% dt of Tchain data
dt=xx2.yday(2)-xx2.yday(1);

figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.09, 0.06, 0.1, 0.09, 2,2);

allws=[0.1 0.25 0.5 0.75];

for whax=1:4 % do for 4 different profiling speeds
    
    clear w_samp Tsamp Tsamp_day nt Ntsamp
    
    % specify vertical sampling velocity
    w_samp=allws(whax)
    
    % time for 1 profile
    Tsamp=range(xx2.z)/w_samp; % sec
    Tsamp_day=Tsamp/86400; % days
    nt=round(Tsamp_day/dt); % ~ # real profiles in 1 resamp profile
    
    Ntsamp=floor(length(xx2.yday)/nt)-1;
    
    axes(ax(whax))
    semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)
    hold on
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    xlabel('<\epsilon>','fontsize',16)
    
    % Repeat a bunch of times, shifting start time by 2 mins to get range of possible results
    eall=nan*ones(length(xx2.z),nt);
    for qt=1:nt
        clear idre
        idre=(1+qt): nt : Ntsamp*nt;
        semilogx(nanmean(xx2.eps(:,idre),2),xx2.z,'.','color',0.75*[1 1 1],'linewidth',1)
        eall(:,qt)=nanmean(xx2.eps(:,idre),2);
    end
    
    semilogx(nanmean(xx2.eps,2),xx2.z,'k','linewidth',2)
    semilogx(nanmean(eall,2),xx2.z,'m--','linewidth',1)
    xlim([1e-9 1e-6])
    ylim([xx2.z(1) xx2.z(end)])
    title(['Tchain ' num2str(whmoor) ' - w=' num2str(w_samp) ' m/s'],'fontsize',16)
    
end

shg

if saveplot==1
    %
    fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['Tchain' num2str(whmoor) '_UnderSamp_4speeds'])
    addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
    %export_fig(fname,'-pdf')
    %
end
%%