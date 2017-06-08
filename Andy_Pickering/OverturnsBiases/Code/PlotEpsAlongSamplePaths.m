function [eps_prof]=PlotEpsAlongSamplePaths(REsamp,xl,whdir,cmap,cvec,whvar)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% General function to plot resampled epsilon along sampling paths.
%
% INPUT
% REsamp : Resampled data structure made w/ ResampleFieldGeneral.m
% xl     : Time range to plot
% whdir  : Direction of sampling to plot ('up' or down)
% cmap   : Colormap to use
% cvec   : Caxis limits
% whvar  : Which variable to plot
%
%
% Modified from Plot_Eps_AlongSamplePaths.m
%
%---------------------
% May 29, 2015 - A. Pickering
% 19 June - Allow other variables to be plotted
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
% plot eps along sampling path for each realization

eps_prof=nan*ones(length(REsamp.z),REsamp.Nshift);

for whc=1:REsamp.Nshift
    
    clear idt ec tc zc thdir dir_corr
    
    % find data in the time range we want
    idt=isin(REsamp.tgrid(whc,:),xl);
    
    % find direction of first profile in realization
    thdir=sign( REsamp.zsamp(2,idt(1),whc)-REsamp.zsamp(1,idt(1),whc));
    
    if strcmp(whdir,'up')
        dir_corr=-1;
    else
        dir_corr=1;
    end
    
    % get only up or down paths (since they will overlap)
    if thdir*dir_corr==1 % 1st profile is down
        ec=squeeze(REsamp.(whvar)(:,idt(1:2:end),whc));
        tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
        zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
    elseif thdir*dir_corr==-1 % 1st profile is up
        ec=squeeze(REsamp.(whvar)(:,idt(2:2:end),whc));
        tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
        zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
    end
    
    % save the mean depth profiles
    eps_prof(:,whc)=nanmean(ec,2);
    
    if strcmp(whdir,'up')
        % REsamp.eps always increasing depth
        ec=flipud(ec);
    end
    
    % plot each point along sampling path
    for ii=1:numel(ec)
        % find closest color
        [val,I]=nanmin(abs(log10(ec(ii))-cvec));
        if I~=1
            plot(tc(ii),zc(ii),'o','color',cmap(I,:),'linewidth',2,'markersize',2)
            hold on
        end
    end
    
end % whc

cb=colorbar('East');cb.AxisLocation='out';
axis ij
%xlim(xl);ylim(yl)
caxis([cvec(1) cvec(end)])

return
%%