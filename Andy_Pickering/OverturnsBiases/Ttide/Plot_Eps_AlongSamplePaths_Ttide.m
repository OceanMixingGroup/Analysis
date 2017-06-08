%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Eps_AlongSamplePaths_Ttide.m
%
% Plot epsilon along actual resampling paths. For mooring `T2' from T-tide
%
%
% Modified from Plot_Eps_AlongSamplePaths.m
%
% 10 Mar. 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

plotiso=1
minOT=10
yl=[1600 1900]
saveplot=0

xl=[19.2 19.45];str='NearDay19' ; cl=[-8 -4.9]; xlprof=[1e-8 1e-6]
%xl=[20.3 20.45];str='NearDay20' ; cl=[-8 -5.5]; xlprof=[1e-8 1e-6]
%xl=[20.8 21.02];str='NearDay21' ; cl=[-8 -5]; xlprof=[1e-8 1e-6]

minOT=10

useOrig=0
useGridded=1 ; dzgrid=10 ; dtgrid=2;

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

% load `true' data
if useOrig==1
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '.mat'])
elseif useGridded==1
    load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/T2_RecomputedEps_MinOT_' num2str(minOT) '_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
else
end
diso=3

for testnum=9%[2 9]
    
    %~~~ load resampled dataset
    clear fname  M N dt tc ei idsamp ep REsamp
    DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Ttide/Data/T2/Test' num2str(testnum)]
    
    if useOrig==1
        fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
    elseif useGridded==1
        fname=fullfile(DatDir,['T2_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases_Gridded_dz_' num2str(dzgrid) 'm_dt_' num2str(dtgrid) 'min']);
    else
    end
    
    load(fname)
    
    figure(1);clf
    agutwocolumn(1)
    wysiwyg
    ax = MySubplot4(0.1, 0.03, 0.04, 0.06, 0.1, 0.03, 2,3)
    
    tm=nanmean(xx2.T,2);
    tm=tm(~isnan(tm));
    
    id=isin(xx2.yday,xl);
    
    % plot true T-chain data
    axes(ax(2))
    ezpc(xx2.yday(id),xx2.z,log10(xx2.eps(:,id)))
    caxis(cl)
    
    % plot isotherms
    if plotiso==1
        hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
    end
    
    hold on
    % plot example resampling path
    plot(REsamp.tsamp(:,:,1),REsamp.zsamp(:,:,1),'w')
 %    Nshift=size(REsamp.eps,3);
%    plot(REsamp.tsamp(:,:,round(Nshift/2)),REsamp.zsamp(:,:,round(Nshift/2)),'w')
        
    cmap=flipud(hot);
    colormap([0.75*[1 1 1] ; cmap])
    %    xtloff
    xlim(xl);ylim(yl)
    title(['T2, yday ' num2str(xl(1)) ' : ' num2str(xl(2)) ' , w=' num2str(REsamp.w_samp) 'm/s' ],'interpreter','none')
    SubplotLetterMW('\epsilon true')
    %
    
    axes(ax(1))
    semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    ylim(yl)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    %
    cvec=linspace(cl(1),cl(2),length(cmap));
    
    Nshift=size(REsamp.eps,3);
    eps_prof_all_up=nan*ones(length(REsamp.z),Nshift);
    eps_prof_all_dn=nan*ones(length(REsamp.z),Nshift);
    
    axes(ax(4))
    
    % plot isotherms
    if plotiso==1
        % hold on
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
        hold on
    end
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
        
        for ii=1:numel(ec)
            [val,I]=nanmin(abs(log10(ec(ii))-cvec));
            if I~=1
                plot(tc(ii),zc(ii),'o','color',cmap(I,:),'linewidth',2,'markersize',2)
                hold on
            end
        end
        
        eps_prof_all_up(:,whc)=nanmean(ec,2);
        
    end
    
    
    cb=colorbar('East');
    cb.AxisLocation='out'
    axis ij
    caxis([cvec(1) cvec(end)])
    xlim(xl)
    xlim(xl);ylim(yl)
    
    axes(ax(3))
    semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    ylim(yl)
    hold on
    semilogx(nanmean(eps_prof_all_up,2),REsamp.z)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    legend('true','sampled','location','best')
    
    %
    axes(ax(6))
    % plot isotherms
    if plotiso==1
        contour(xx2.yday(id),xx2.z,xx2.T(:,id),tm(1:diso:end),'color',0.1*[1 1 1])
        hold on
    end
    
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
        
        for ii=1:numel(ec)
            [val,I]=nanmin(abs(log10(ec(ii))-cvec));
            if I~=1
                plot(tc(ii),zc(ii),'o','color',cmap(I,:),'linewidth',2,'markersize',2)
                hold on
            end
        end
        
        eps_prof_all_dn(:,whc)=nanmean(flipud(ec),2);
        
    end
    
    
    caxis([cvec(1) cvec(end)])
    axis ij
    xlim(xl);ylim(yl)
    
    axes(ax(5))
    semilogx(nanmean(xx2.eps(:,id),2),xx2.z,'k','linewidth',2)
    xlim(xlprof)
    hold on
    semilogx(nanmean(eps_prof_all_dn,2),REsamp.z)
    axis ij
    grid on
    ylabel('Depth','fontsize',16)
    xlabel('<\epsilon>','fontsize',16)
    ylim(yl)
    
    linkaxes(ax,'y')
    linkaxes(ax([2 4 6]))
    
    %
    if saveplot==1
        figdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases'
        addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
        fname=fullfile(figdir,['Ttide_T2_AlongPath_' str '_Test_' num2str(testnum) ])
        export_fig(fname,'-pdf','-r200')
    end
    
    
end % testnum
%%