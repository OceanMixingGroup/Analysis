%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotDepthLoops.m
%
% For EQ14 CTD-chipod data.
%
%-------------------
% 04/25/16 - A.Pickering
% 06/28/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/

Params.extra_z=2;
Params.wthresh=0.4;

datadir='/Users/Andy/Cruises_Research/ChiPod/EQ14/Data/Chipod/proc/SN1001/cal/'
Flist=dir(fullfile(datadir,'*.mat'))

pball=nan*ones(1,length(Flist));
pbchi=pball;
pgh=pball;
%
for icast=62%:length(Flist)
    
    try
        
        load(fullfile(datadir,Flist(icast).name))
        %
        iz=find(C.ctd.raw.p>20);
        praw=C.ctd.raw.p(iz);
        dtraw=C.ctd.raw.datenum(iz);
        
        data2.p=praw;
        data2.datenum=dtraw;
        
        % remove loops in CTD data
        clear datau2 bad_inds tmp
        [datau2,bad_inds] = ctd_rmdepthloops(data2,Params.extra_z,Params.wthresh);
        tmp=ones(size(datau2.p));
        tmp(bad_inds)=0;
        
        %
        figure(1);clf
        
        ax1=subplot(211);
        plot(dtraw,praw)
        hold on
        plot(dtraw(bad_inds),praw(bad_inds),'r.')
        axis ij
        datetick('x',15)
        shg
        
        ax2=subplot(212);
        plot(C.datenum,C.T1P,'b')
        ylim(1.5*[-1 1])
        datetick('x',15)
        
        linkaxes([ax1 ax2],'x')
        
        C.is_good_data=interp1(datau2.datenum,tmp,C.datenum,'nearest');
        clear ib_loop Nloop
        ib_loop=find(C.is_good_data==0);
        Nloop=length(ib_loop);
        disp(['\n  ' num2str(round(Nloop/length(C.datenum)*100)) ' percent of points removed for depth loops ']);
        
        % percent of ctd points that are bad
        pball(icast)=length(unique(bad_inds))/length(praw)*100;
        
        % percent of chipod data thrown out
        pbchi(icast)=Nloop/length(C.datenum)*100;
        
        %
        C.T1Praw=C.T1P;
        C.T1P(ib_loop)=nan;
        
        
         axes(ax2)
         hold on
        plot(C.datenum,C.T1P,'r')
        
        clear todo_inds Nwindows
        Params.nfft=128;
        [todo_inds,Nwindows]=MakeCtdChiWindows(C.T1P,Params.nfft);
        
        %~ make 'avg' structure for the processed data
        clear avg
        avg=struct();
        avg.Params=Params;
        tfields={'datenum','P','N2','dTdz','fspd','T','S','P','theta','sigma',...
            'chi1','eps1','KT1','TP1var'};
        for n=1:length(tfields)
            avg.(tfields{n})=NaN*ones(Nwindows,1);
        end
        
        avg.samplerate=1./nanmedian(diff(C.datenum))/24/3600;
        
        % get average time, pressure, and fallspeed in each window
        for iwind=1:Nwindows
            clear inds
            inds=todo_inds(iwind,1) : todo_inds(iwind,2);
            avg.datenum(iwind)=nanmean(C.datenum(inds));
            avg.P(iwind)=nanmean(C.P(inds));
            avg.fspd(iwind)=nanmean(C.fspd(inds));
        end
        
        %~~ plot histogram of avg.P to see how many good windows we have in
        % each 10m bin
        figure(2);clf
        hi=histogram(avg.P,0:10:nanmax(avg.P),'Edgecolor','none');
        hi.Orientation='Horizontal';axis ij;
        ylabel('P [db]')
        xlabel('# good data windows')
        %title([whSN ' cast ' castStr ' - ' C.castdir 'cast'],'interpreter','none')
        %
        
        ig=find(hi.Values>5);
        pgh(icast)=length(ig)/length(hi.Values);
        
    end % try
    
end % icast
%%
cast_id=1:length(Flist);

%%

figure(1);clf
agutwocolumn(0.5)
wysiwyg

subplot(131)
histogram(pball)

subplot(132)
histogram(pbchi)

subplot(133)
histogram(pgh*100)

xlim([0 100])

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(311)
plot(cast_id,pball,'o-')
grid on
ylim([0 50])
title('% CTD data bad')

subplot(312)
plot(cast_id,pbchi,'o-')
grid on
ylim([0 50])
title('% chi data bad')

subplot(313)
plot(cast_id,pgh*100,'o-')
xlabel('Cast ID')
ylabel('% Good 10m bins')
grid on

%%

figure(1);clf

ax1=subplot(121)
plot(dtraw,praw,'k')
hold on
plot(dtraw(bad_inds),praw(bad_inds),'r.')
axis ij
datetick('x','keeplimits')
xlim([datenum(2014,12,2,5,50,41) datenum('02-Dec-2014 05:51:51') ])
datetick('x','keeplimits')
ylim([280 330])
grid on
ylabel('Pres. [db]','fontsize',16)

ax2=subplot(122)
hi=histogram(avg.P,0:10:nanmax(avg.P))%,'Edgecolor','none');
hi.Orientation='Horizontal';axis ij;
ylabel('P [db]')
xlabel('# Good windows per Bin')
xlim([0 20])

linkaxes([ax1 ax2],'y')

%%

SetPaperFigPath
figname='DepthLoopExample'
print(fullfile(figdir,figname),'-dpng')

%%

fcut=86400/3;
p_low=MyLowpass(dtraw,praw,2,fcut);

figure(1);clf

ax1=subplot(211)
plot(dtraw,sw_dpth(p_low,0),'k')
hold on
%plot(dtraw,p_low,'g')
hold on
%hloop=plot(dtraw(bad_inds),praw(bad_inds),'r.')
hloop=plot(dtraw(bad_inds),sw_dpth(p_low(bad_inds),0),'r.')
axis ij
datetick('x','keeplimits')
xlim([datenum(2014,12,2,5,52,00) datenum('02-Dec-2014 05:52:39') ])
datetick('x','keeplimits')
grid on
%ylabel('Pres. [db]','fontsize',16)
ylabel('Depth [m]','fontsize',16)
legend(hloop,'Depth loops','location','best')

ax2=subplot(212)
hf=plot(C.datenum,C.T1Praw,'color','r') %0.5*[1 1 1])
hold on
plot(C.datenum,C.T1P,'k')
ylim(0.6*[-1 1])
%legend('raw','filtered')
legend(hf,'flagged')
datetick('x','keeplimits')
xlim([datenum(2014,12,2,5,52,00) datenum('02-Dec-2014 05:52:39') ])
datetick('x','keeplimits')
ylabel('dT/dt [^oCs^{-1}]','fontsize',16)
grid on
xlabel(['Time '])

linkaxes([ax1 ax2],'x')

%%
SetPaperFigPath
figname='DepthLoopExample2TS'
print(fullfile(figdir,figname),'-dpng')

%%
