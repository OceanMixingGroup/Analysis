%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Eps_AlongSamplePaths.m
%
% Plot epsilon along actual resampling paths.
%
%
% Used to be part of ExamineShortSectionsResamp.m
%
% 9 Mar. 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 5 Mar. - plot along actual sampling path


clear ; close all

plotiso=1
%whmoor=1 ; yl=[500 1500]
%whmoor=2 ; yl=[500 1500]
whmoor=3 ; yl=[800 2100]
%whmoor=4 ; yl=[900 2200]
minOT=50

saveplot=1

% tchain 3
%xl=[168.25 168.7];str='NearDay168' ; cl=[-8 -4.2]; xlprof=[1e-8 1e-5]
%xl=[168.3440  168.4903]

% tchain 4
%xl=[169.3 169.6];str='NearDay169' ; cl=[-8 -4.5]; xlprof=[1e-8 1e-5]

% tchain3
xl=[169.35 169.55];str='NearDay169' ; cl=[-8 -4]; xlprof=[1e-8 1e-5]

% tchain1
%xl=[181.3 181.6]; str='NearDay181' ; cl=[-8 -4.7]; xlprof=[1e-9 1e-6]

% tchain2
%xl=[181.2 181.5]; str='NearDay181' ; cl=[-8 -4]; xlprof=[1e-9 1e-5]

% tchain2
%xl=[180.2 180.4]; str='NearDay180' ; cl=[-8 -4]; xlprof=[1e-9 1e-5]

% tchain4
%xl=[181.2 181.55]; str='NearDay181' ; cl=[-8 -4.2]; xlprof=[1e-9 1e-5]

% tchain3
%xl=[182.2 182.5]; str='NearDay182' ; cl=[-8 -3.6]; xlprof=[1e-8 1e-4]

% tchain4
%xl=[182.1 182.6]; str='NearDay182' ; cl=[-8 -4]; xlprof=[1e-9 1e-5]

% tchain3
%xl=[183.3 183.5] ; str='NearDay183' ; cl=[-8 -3.6]; xlprof=[1e-8 1e-4]

% tchain4
%xl=[183.2 183.5]; str='NearDay183' ; cl=[-8 -4.4]; xlprof=[1e-8 1e-5]

% tchain3
%xl=[184.2 184.7]; str='NearDay184' ; cl=[-8 -4]; xlprof=[1e-8 1e-4]

%xl=[185.2 185.7]; str='NearDay185' ; cl=[-8 -4]; xlprof=[1e-8 1e-4]

%xl=[197.3 197.5];str='NearDay197' ;  cl=[-8 -3.8]; xlprof=[1e-8 1e-5]

%xl=[198.3 198.6] ; str='NearDay198' ; cl=[-8 -4.5]; xlprof=[1e-8 1e-5]

% Tchain4
%xl=[198.25 198.55] ; str='NearDay198' ; cl=[-8 -4.5]; xlprof=[1e-8 1e-6]

%xl=[210.25 210.4];str='NearDay210' ; cl=[-8 -3.8]; xlprof=[1e-8 1e-4]

%xl=[211.2 211.45];str='NearDay211' ; cl=[-8 -3.5]; xlprof=[1e-8 1e-4]

%

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

diso=15

for testnum=3%[2 3]
    
    %~~~ load resampled dataset
    clear fname  M N dt tc ei idsamp ep REsamp
    
    load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
    
    figure(1);clf
    agutwocolumn(1)
    wysiwyg
    %ax = MySubplot(0.1, 0.03, 0.1, 0.06, 0.1, 0.01, 1,3);
    ax = MySubplot4(0.1, 0.03, 0.04, 0.06, 0.1, 0.03, 2,3)
    
    tm=nanmean(xx2.T,2);
    tm=tm(~isnan(tm));
    
    id=isin(xx2.yday,xl);
    
    whvar='eps' ; cl=[-8 -4]
%     whvar='Otnsq_out' ; %cl=[]
%     whvar='n2' ; cl=[-7.5 -4.5]
%     whvar='Lot' ; cl=[1.5 3]
    

%~~ plot mean depth profile of true epsilon
    axes(ax(1))
    semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    ylim(yl)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
     SubplotLetterMW('a)')
    
    %~~ plot true T-chain data
    axes(ax(2))
    ezpc(xx2.yday(id),xx2.z,real(log10(xx2.(whvar)(:,id))))
    caxis(cl)
    hold on
    
    % plot isotherms
    if plotiso==1
        hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
    end
    
    % plot example resampling path
    plot(REsamp.tsamp(:,:,1),REsamp.zsamp(:,:,1),'w')
    
    % make colormap for epsilon
    cmap=flipud(hot);
    colormap([0.75*[1 1 1] ; cmap])
    
    xlim(xl);ylim(yl)
    ytloff
    title(['Tchain ' num2str(whmoor) ', yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ' , w=' num2str(REsamp.w_samp) 'm/s' ],'interpreter','none')
    SubplotLetterMW('d) \epsilon true')
        
    
    % color vector for plotting the resampled epsilon along paths
    cvec=linspace(cl(1),cl(2),length(cmap));
    
    Nshift=size(REsamp.eps,3);
    eps_prof_all_up=nan*ones(length(REsamp.z),Nshift);
    eps_prof_all_dn=nan*ones(length(REsamp.z),Nshift);
    
    %~~ plot DOWNward sampling paths
    axes(ax(4))
    
    % plot isotherms
    if plotiso==1
        hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
    end
    
    %     % plot eps along sampling path for each realization
    %     for whc=1:Nshift
    %         clear idt ec tc zc whdir
    %         % find data in the time range we want
    %         idt=isin(REsamp.tgrid(whc,:),xl);
    %         % find direction of first profile in realization
    %         whdir=sign( REsamp.zsamp(2,idt(1),whc)-REsamp.zsamp(1,idt(1),whc));
    %
    %         % get only up or down paths (since they will overlap)
    %         if whdir==1 % 1st profile is down
    %             ec=squeeze(REsamp.eps(:,idt(1:2:end),whc));
    %             tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
    %             zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
    %         elseif whdir==-1 % 1st profile is up
    %             ec=squeeze(REsamp.eps(:,idt(2:2:end),whc));
    %             tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
    %             zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
    %         end
    %
    %         for ii=1:numel(ec)
    %             [val,I]=nanmin(abs(log10(ec(ii))-cvec));
    %             if I~=1
    %                 plot(tc(ii),zc(ii),'o','color',cmap(I,:),'linewidth',2,'markersize',2)
    %                 hold on
    %             end
    %         end
    %
    %         eps_prof_all_dn(:,whc)=nanmean(ec,2);
    %
    %     end
    %
    %     cb=colorbar('East');cb.AxisLocation='out';
    %     axis ij
    
    whdir='down'
    eps_prof_all_dn=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec,whvar);
    
    xlim(xl);ylim(yl)
    ytloff
    caxis([cvec(1) cvec(end)])
        SubplotLetterMW('e)')
        
    % plot mean profile for downward sampling paths
    axes(ax(3))
    semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    ylim(yl)
    hold on
    semilogx(nanmean(eps_prof_all_dn,2),REsamp.z,'linewidth',2)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    legend('true','sampled','location','best')
    SubplotLetterMW('b)')
    
    %~ plot upward sampling paths
    axes(ax(6))
    
    % plot isotherms
    if plotiso==1
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
        hold on
    end
    
    %     for whc=1:Nshift
    %         clear idt ec tc zc whdir
    %         idt=isin(REsamp.tgrid(whc,:),xl);
    %
    %         whdir=sign( REsamp.zsamp(2,idt(1),whc)-REsamp.zsamp(1,idt(1),whc));
    %
    %         % plot only up or down (since they will overlap)
    %         if whdir==-1 % 1st profile is up
    %             ec=squeeze(REsamp.eps(:,idt(1:2:end),whc));
    %             tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
    %             zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
    %         elseif whdir==1 % 1st profile is down
    %             ec=squeeze(REsamp.eps(:,idt(2:2:end),whc));
    %             tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
    %             zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
    %         end
    %
    %         % REsamp.eps always increasing depth
    %         ec=flipud(ec);
    %
    %         % plot each point along path
    %         for ii=1:numel(ec)
    %             % find closest color
    %             [val,I]=nanmin(abs(log10(ec(ii))-cvec));
    %             if I~=1
    %                 plot(tc(ii),zc(ii),'o','color',cmap(I,:),'linewidth',2,'markersize',2)
    %                 hold on
    %             end
    %         end
    %
    %         eps_prof_all_up(:,whc)=nanmean(flipud(ec),2);
    %
    %     end % whc
    %
    whdir='up'
    eps_prof_all_up=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec,whvar);    
    caxis([cvec(1) cvec(end)])
    xlim(xl);ylim(yl)
    ytloff
    xlabel('Yearday 2011','fontsize',16)
        SubplotLetterMW('f)')
        
    % plot mean depth profile for upward sampling paths
    axes(ax(5))
    semilogx(nanmean(xx2.(whvar)(:,id),2),xx2.z,'k','linewidth',2)
    hold on
    semilogx(nanmean(eps_prof_all_up,2),REsamp.z,'linewidth',2)
    xlim(xlprof);ylim(yl)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    xlabel('<\epsilon>','fontsize',16)
        SubplotLetterMW('c)')
    linkaxes(ax,'y')
    linkaxes(ax([2 4 6]))
    
    %
    if saveplot==1
        figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
        addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
        fname=fullfile(figdir,['Tchain' num2str(whmoor) '_AlongPath_' str '_Test_' num2str(testnum) ])
        export_fig(fname,'-pdf','-r200')
    end
    
    
end % testnum
%
%
%
%
%% Plot just depth profiles

xlprof=[1e-9 1e-4]
xl=[180 200]

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
load( fullfile( 'Data' ,['Tchain' num2str(whmoor)], ['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]) )

diso=15

Nshift=REsamp.Nshift

for testnum=3%[2 3]
    %~~~ load resampled dataset
    clear fname  M N dt tc ei idsamp ep REsamp
    
    load (fullfile(cd,'Data',['Tchain' num2str(whmoor)],['Test' num2str(testnum)],['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases']))
    
    
    figure(1);clf
    agutwocolumn(0.7)
    wysiwyg
    %    ax = MySubplot4(0.1, 0.03, 0.04, 0.06, 0.1, 0.03, 1,2)
    
    id=isin(xx2.yday,xl);
    
    ht=    semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    ylim(yl)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    %
    clear eps_prof_dn
    for whc=1:Nshift
        clear idt ec tc zc whdir
        idt=isin(REsamp.tgrid(whc,:),xl);
        
        whdir=sign( REsamp.zsamp(2,idt(1),whc)-REsamp.zsamp(1,idt(1),whc));
        
        % plot only up or down (since they will overlap)
        if whdir==1
            ec=squeeze(REsamp.eps(:,idt(1:2:end),whc));
            tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
            zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
        elseif whdir==-1
            ec=squeeze(REsamp.eps(:,idt(2:2:end),whc));
            tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
            zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
        end
        
        eps_prof_all_dn(:,whc)=nanmean(ec,2);
        
    end
    
    hold on
    hd=    semilogx(nanmean(eps_prof_all_dn,2),REsamp.z)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    %    legend('true','sampled','location','best')
    
    clear eps_prof_up
    for whc=1:Nshift
        clear idt ec tc zc whdir
        idt=isin(REsamp.tgrid(whc,:),xl);
        
        whdir=sign( REsamp.zsamp(2,idt(1),whc)-REsamp.zsamp(1,idt(1),whc));
        
        % plot only up or down (since they will overlap)
        if whdir==-1
            ec=squeeze(REsamp.eps(:,idt(1:2:end),whc));
            tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
            zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
        elseif whdir==1
            ec=squeeze(REsamp.eps(:,idt(2:2:end),whc));
            tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
            zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
        end
        
        % REsamp.eps always increasing depth
        ec=flipud(ec);
        eps_prof_all_up(:,whc)=nanmean(flipud(ec),2);
        
    end
    
    xlim(xlprof)
    hold on
    hu=    semilogx(nanmean(eps_prof_all_up,2),REsamp.z)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    xlabel('<\epsilon>','fontsize',16)
    ylim(yl)
    
    legend([ht hu hd],'True','up','down')
    %
    if saveplot==1
        figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
        addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
        %        fname=fullfile(figdir,['Tchain' num2str(whmoor) '_AlongPath_' str '_Test_' num2str(testnum) ])
        %       export_fig(fname,'-pdf','-r200')
    end
    
    
end % testnum
%%
