%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Some misc. code trying to come up with a test to find biases in overturns
% by computing vertical velocity from isotherms and finding where it is in
% same direction as profiler.
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

w_est=xx2.Lttot.*sqrt(xx2.n2);
w_est=real(w_est);

% use mean N profile
n2m=nanmean(xx2.n2,2);
w_est2=xx2.Lttot.*repmat(sqrt(n2m),1,length(xx2.yday));

figure(1);clf
histogram(w_est(:),'Normalization','pdf')
hold on
histogram(w_est2(:),'Normalization','pdf')
xlim([0 0.5])


%% 23 November 2015


% Use Lt*N, which is a vertical velocity scale for the overturn?

clear ; close all

whmoor=3; minOT=50; 
testnum=3 % 2=0.15
%testnum=2 % w=0.5
%testnum=1 % w=0.25
%testnum=4 % w=0.75


cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true data'
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
%


wt_true=xx2.Lt_each.*sqrt(xx2.Otnsq_each);

eps1=[]
eps2=[]
wt2=[];
for whc=1:REsamp.Nshift
    clear wt_samp ib
    wt_samp=REsamp.Lttot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
    wt_samp=real(wt_samp);
    wt2=[wt2 ; wt_samp(:)];
    
    clear ec
    ec=squeeze(REsamp.eps(:,:,whc));
    eps1=[eps1 nanmean(ec,2)];
    ib=find(wt_samp> (0.75*REsamp.w_samp)  );
    ec(ib)=nan;
    eps2=[eps2 nanmean(ec,2)];
    
end
%
Nm='pdf'
Ds='stair'
figure(1);clf
ht=histogram(wt_true,'Normalization',Nm,'DisplayStyle',Ds)
hold on
%hs=histogram(wt_samp,'BinEdges',ht.BinEdges,'Normalization',Nm,'DisplayStyle',Ds)
hs=histogram(wt2,'BinEdges',ht.BinEdges,'Normalization',Nm,'DisplayStyle',Ds)
legend('true','samp')
freqline(REsamp.w_samp)
ylabel(Nm)
xlabel('w_{turb}')
grid on

% Nan out resampled epsilon where w_t>w_samp, and plot mean profile

figure(2);clf
semilogx(nanmean(xx2.eps,2),xx2.z)
hold on
semilogx(nanmean(eps1,2),REsamp.z)
semilogx(nanmean(eps2,2),REsamp.z)
axis ij 
grid on
legend('true','samp','samp corr')
title(['Tchain ' num2str(whmoor) ' - W_{samp}=' num2str(REsamp.w_samp)])

%% May 29, 2015 - try new method

% Use Lt*N, which is a vertical velocity scale for the overturn?

clear ; close all

plotiso=1
whmoor=3
minOT=50;
testnum=3

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true data'
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))


for whc=1:REsamp.Nshift
    
    % compute Lt*N and Nan data where this exceeds sampling speed
    clear ww
    ww=REsamp.Lttot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
    ww=real(ww);
    
    clear idb eps1 eps2
    idb=find(ww>(1*REsamp.w_samp));
    eps1=REsamp.eps(:,:,whc);
    eps2=eps1;
    eps2(idb)=nan;
    
    %
    figure(3);clf
    semilogx(nanmean(eps1,2),REsamp.z)
    hold on
    semilogx(nanmean(eps2,2),REsamp.z)
    semilogx(nanmean(xx2.eps,2),xx2.z)
    xlim([1e-9 1e-5])
    legend('samp','corr','true','location','best')
    axis ij
    
    pause(0.5)
    
end

%% see if we are actually getting the bad times or just getting rid of big overturns..

clear ; close all

plotiso=1
whmoor=3
minOT=50
testnum=1

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load 'true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))


%for whc=1:REsamp.Nshift
%xl=[169.2 169.7]
xl=[169.35 169.55]

whc=1
clear ww
ww=REsamp.Lttot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
ww=real(ww);
% figure(1);clf
% ezpc(REsamp.tgrid(whc,:),REsamp.z,ww)
% colorbar
% caxis([0 REsamp.w_samp])

clear idb eps1 eps2
idb=find(ww>(0.75*REsamp.w_samp));
eps1=REsamp.eps(:,:,whc);
eps2=eps1;
eps2(idb)=nan;
%
figure(2);clf
ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.03, 1,2);

axes(ax(1))
ezpc(REsamp.tgrid(whc,:),REsamp.z,log10(eps1))
colorbar
caxis([-10 -5])
xlim(xl)
hold on

axes(ax(2))
ezpc(REsamp.tgrid(whc,:),REsamp.z,log10(eps2))
colorbar
caxis([-10 -5])
xlim(xl)

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])

%
% linkaxes(ax)
%
% figure(3);clf
% semilogx(nanmean(eps1,2),REsamp.z)
% hold on
% semilogx(nanmean(eps2,2),REsamp.z)
% semilogx(nanmean(xx2.eps,2),xx2.z)
% xlim([1e-9 1e-5])
% legend('samp','corr','true','location','best')
% axis ij

%pause(0.5)

%end

%% Another Different version - 17 Feb
% contour average ensemble

clear ; close all

plotiso=1
whmoor=3
minOT=50
testnum=1

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%load( fullfile( 'Data' , ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

%~~ load resampled dataset
clear fname x_resamp
%load (fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))

REsamp.t=REsamp.data_resamp;
%%
% xl=[168 170]
% xl=[168.2 168.7]
xl=[169.2 169.7]
% xl=[170.2 170.7]
%xl=[178.2 178.5]
%xl=[179 179.6]
%xl=[179.6 180.2]
%xl=[181 181.7]
%xl=[181 182.2]
%xl=[182.2 182.7]
%xl=[180 185]

%xl=[183.25 183.65]

%xl=[183.3 183.55]
%xl=[183.27 183.5]
%xl=[184.2 184.7]
%xl=[210 212]
%xl=[210 210.6]
%xl=[210.25 210.4]
%xl=[195 200]
%xl=[205 215]

%
whc=1

diso=15
cl=[-9 -3]

figure(1);clf
agutwocolumn(0.7)
wysiwyg


tm=nanmean(xx2.T,2);
tm=tm(~isnan(tm));

id=isin(xx2.yday,xl);

% plot true T-chain data
ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
cb.FontSize=14;
caxis(cl)
SubplotLetterMW('\epsilon true')

% plot isotherms
if plotiso==1
    hold on
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'k')
end

cmap=flipud(hot);
colormap([0.75*[1 1 1] ; cmap])
ylabel('Depth','fontsize',16)

% contour resampled isotherms
idsamp=isin(REsamp.tgrid(whc,:),xl);
contour(REsamp.tgrid(whc,idsamp),REsamp.z,REsamp.t(:,idsamp,whc),tm(1:diso:end),'w')
plot(REsamp.tsamp(:,idsamp,whc),REsamp.zsamp(:,idsamp,whc),'b')

title(['Tchain' num2str(whmoor) ' w=' num2str(REsamp.w_samp) 'm/s'])

shg

%%

clc
whp=2

for whc=1:5:size(REsamp.eps,3)
    
    %whc=52
    
    clear id tm zm Nprof  idsamp
    
    id=isin(xx2.yday,xl);
    idsamp=isin(REsamp.tgrid(whc,:),xl);
    
    figure(1);clf
    ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
    hold on
    cb=colorbar;
    cb.Label.String='log_{10}\epsilon';
    cb.FontSize=14;
    caxis(cl)
    SubplotLetterMW('\epsilon true')
    
    tm=nanmean(xx2.T,2);
    zm=xx2.z(~isnan(tm));
    tm=tm(~isnan(tm));
    % plot isotherms
    contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'k')
    hline(zm(1:diso:end),'--')
    cmap=flipud(hot);
    colormap([0.75*[1 1 1] ; cmap])
    ylabel('Depth','fontsize',16)
    %
    
    Nprof=length(idsamp)
    
    % plot the sampling paths we are working with
    plot(REsamp.tsamp(:,idsamp(whp:whp+1),whc),REsamp.zsamp(:,idsamp(whp:whp+1),whc),'w')
    %    mindiff=0.1;
    dzdt=nan*ones(length(tm),1);
    contour(REsamp.tgrid(whc,idsamp),REsamp.z,REsamp.t(:,idsamp,whc),tm(1:diso:end),'w')
    title(['Tchain' num2str(whmoor) ' w=' num2str(REsamp.w_samp) 'm/s'])
    
    % get temperature profiles for these two
    clear prof1 prof2
    prof1=REsamp.t(:,idsamp(whp),whc);
    prof2=REsamp.t(:,idsamp(whp+1),whc);
    
    % mean temperature profile
    clear tm
    tm=nanmean(xx2.T,2);
    
    % find instantaneous depths of each temp from the mean profile
    for wht=1:length(tm)
        clear I_1 I_2 I_11 I_22  dz dt
        [val_1,I_1]=nanmin(abs(prof1-tm(wht)));
        [val_2,I_2]=nanmin(abs(prof2-tm(wht)));
        
        [val,I_11]=nanmin(abs(REsamp.z(I_1) - REsamp.zsamp(:,idsamp(whp),whc)  ) );
        [val,I_22]=nanmin(abs(REsamp.z(I_2) - REsamp.zsamp(:,idsamp(whp+1),whc)  ) );
        
        % plot the depths
        plot(REsamp.tsamp(I_11,idsamp(whp),whc),REsamp.zsamp(I_11,idsamp(whp),whc),'o')
        plot(REsamp.tsamp(I_22,idsamp(whp+1),whc),REsamp.zsamp(I_22,idsamp(whp+1),whc),'o')
        
        % draw the line connecting them
        line([REsamp.tsamp(I_11,idsamp(whp),whc) REsamp.tsamp(I_22,idsamp(whp+1),whc)],[REsamp.zsamp(I_11,idsamp(whp),whc) REsamp.zsamp(I_22,idsamp(whp+1),whc)]);
        
        % compute dz/dt (estimate of w_iso) from these points
        dz=REsamp.zsamp(I_11,idsamp(whp),whc)-REsamp.zsamp(I_22,idsamp(whp+1),whc);
        dt=(REsamp.tsamp(I_22,idsamp(whp+1),whc) - REsamp.tsamp(I_11,idsamp(whp),whc) )*86400;
        dzdt(wht)=dz/dt;
        
        %       pause
        
    end
    %
    figure(2);clf
    plot(dzdt,REsamp.z)
    axis ij
    %freqline(REsamp.w_samp)
    xlim(0.2*[-1 1])
    gridxy
    pause(0.25)
    
end
%
%
%
%% in this version, use profiles before and after current profile so we can compare to the profiler's speed

whp=3
clear w_prof
if nanmean(diff(REsamp.zsamp(:,idsamp(whp),whc)))<0
    w_prof=REsamp.w_samp; % up profile
else
    w_prof=-REsamp.w_samp; % down profile
end

for whc=80%:5:size(REsamp.eps,3)
    
    %whc=52
    
    clear id tm zm Nprof  idsamp
    
    id=isin(xx2.yday,xl);
    idsamp=isin(REsamp.tgrid(whc,:),xl);
    
    figure(1);clf
    
    ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
    cb=colorbar;
    cb.Label.String='log_{10}\epsilon';
    cb.FontSize=14;
    caxis(cl)
    SubplotLetterMW('\epsilon true')
    
    tm=nanmean(xx2.T,2);
    zm=xx2.z(~isnan(tm));
    tm=tm(~isnan(tm));
    % plot isotherms
    if plotiso==1
        hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'k')
        hline(zm(1:diso:end),'--')
    end
    
    cmap=flipud(hot);
    colormap([0.75*[1 1 1] ; cmap])
    ylabel('Depth','fontsize',16)
    %
    
    ip1=idsamp(whp-1);
    ip2=idsamp(whp+1);
    
    Nprof=length(idsamp)
    
    
    plot(REsamp.tsamp(:,ip1:ip2,whc),REsamp.zsamp(:,ip1:ip2,whc),'w')
    %
    mindiff=0.1;
    dzdt=nan*ones(length(tm),1);
    contour(REsamp.tgrid(whc,idsamp),REsamp.z,REsamp.t(:,idsamp,whc),tm(1:diso:end),'w')
    %
    title(['Tchain' num2str(whmoor) ' w=' num2str(REsamp.w_samp) 'm/s'])
    
    clear prof1 prof2
    prof1=REsamp.t(:,ip1,whc);
    prof2=REsamp.t(:,ip2,whc);
    
    tm=nanmean(xx2.T,2);
    for wht=1:length(tm)
        clear I_1 I_2 I_11 I_22  dz dt
        [val_1,I_1]=nanmin(abs(prof1-tm(wht)));
        [val_2,I_2]=nanmin(abs(prof2-tm(wht)));
        
        [val,I_11]=nanmin(abs(REsamp.z(I_1) - REsamp.zsamp(:,ip1,whc)  ) );
        [val,I_22]=nanmin(abs(REsamp.z(I_2) - REsamp.zsamp(:,ip2,whc)  ) );
        
        plot(REsamp.tsamp(I_11,ip1,whc),REsamp.zsamp(I_11,ip1,whc),'o')
        plot(REsamp.tsamp(I_22,ip2,whc),REsamp.zsamp(I_22,ip2,whc),'o')
        line([REsamp.tsamp(I_11,ip1,whc) REsamp.tsamp(I_22,ip2,whc)],[REsamp.zsamp(I_11,ip1,whc) REsamp.zsamp(I_22,ip2,whc)])
        
        dz=REsamp.zsamp(I_11,ip1,whc)-REsamp.zsamp(I_22,ip2,whc);
        dt=(REsamp.tsamp(I_22,ip2,whc) - REsamp.tsamp(I_11,ip1,whc) )*86400;
        dzdt(wht)=dz/dt;
        
        %       pause
        
    end
    %
    figure(2);clf
    plot(dzdt,REsamp.z)
    hold on
    w_prof=sign(nanmean(diff(REsamp.zsamp(:,idsamp(whp),whc))))*REsamp.w_samp;
    freqline(w_prof)
    %    plot(REsamp.w_samp
    axis ij
    %freqline(REsamp.w_samp)
    xlim(0.2*[-1 1])
    gridxy
    pause(0.5)
    
    
    %
    figure(3);clf
    plot(dzdt/w_prof,REsamp.z)
    title('w_{iso}/w_{prof}')
    axis ij
    grid on
    gridxy
    
end

%% another version: calculate displacement from mean depth of each ispycnal



for whc=1:5:size(REsamp.eps,3)
    whp=4
    %whc=50
    clear id tm zm Nprof idsamp
    
    id=isin(xx2.yday,xl);
    idsamp=isin(REsamp.tgrid(whc,:),xl);
    
    figure(1);clf
    
    ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
    cb=colorbar;
    cb.Label.String='log_{10}\epsilon';
    cb.FontSize=14;
    caxis(cl)
    SubplotLetterMW('\epsilon true')
    
    tm=nanmean(xx2.T,2);
    zm=xx2.z(~isnan(tm));
    tm=tm(~isnan(tm));
    % plot isotherms
    dd=10
    if plotiso==1
        hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:dd:end),'k')
        hline(zm(1:dd:end),'--')
    end
    
    cmap=flipud(hot);
    colormap([0.75*[1 1 1] ; cmap])
    ylabel('Depth','fontsize',16)
    %
    %
    ip1=idsamp(whp);
    %ip2=idsamp(whp+1);
    
    Nprof=length(idsamp);
    %
    plot(REsamp.tsamp(:,ip1,whc),REsamp.zsamp(:,ip1,whc),'w')
    %
    title(['Tchain' num2str(whmoor) ' w=' num2str(REsamp.w_samp) 'm/s'])
    
    clear prof1 prof2
    prof1=REsamp.t(:,ip1,whc);
    %prof2=REsamp.t(:,ip2,whc);
    
    %    tm=nanmean(xx2.T,2);
    tm=nanmean(xx2.T(:,id),2);
    diso=nan*ones(length(tm),1);
    for wht=1:length(tm)
        clear I_1 I_2 I_11 I_22  dz dt
        [val_1,I_1]=nanmin(abs(prof1-tm(wht)));
        %    [val_2,I_2]=nanmin(abs(prof2-tm(wht)));
        
        [val,I_11]=nanmin(abs(REsamp.z(I_1) - REsamp.zsamp(:,ip1,whc)  ) );
        %   [val,I_22]=nanmin(abs(REsamp.z(I_2) - REsamp.zsamp(:,ip2,whc)  ) );
        
        %    plot(REsamp.tsamp(I_11,ip1,whc),REsamp.zsamp(I_11,ip1,whc),'o')
        %    plot(REsamp.tsamp(I_22,ip2,whc),REsamp.zsamp(I_22,ip2,whc),'o')
        %    line([REsamp.tsamp(I_11,ip1,whc) REsamp.tsamp(I_22,ip2,whc)],[REsamp.zsamp(I_11,ip1,whc) REsamp.zsamp(I_22,ip2,whc)])
        
        diso(wht)=REsamp.zsamp(I_11,ip1,whc)-REsamp.z(wht);
        %   dt=(REsamp.tsamp(I_22,ip2,whc) - REsamp.tsamp(I_11,ip1,whc) )*86400;
        %   dzdt(wht)=dz/dt;
        
        %       pause
        
    end
    %
    %
    figure(2);clf
    plot(diso,REsamp.z)
    axis ij
    xlim(250*[-1 1])
    gridxy
    
    pause(0.5)
end