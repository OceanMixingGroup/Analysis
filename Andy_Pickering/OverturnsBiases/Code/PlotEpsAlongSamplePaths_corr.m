function [eps_prof]=PlotEpsAlongSamplePaths_corr(REsamp,xl,whdir,cmap,cvec,thresh)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% [eps_prof]=PlotEpsAlongSamplePaths_corr(REsamp,xl,whdir,cmap,cvec)
%
% Function to plot resampled epsilon along sampling paths, with correction
% for biased data based on inferred turbulent velocity scale
%
% Modified from Plot_Eps_AlongSamplePaths.m
%
% See also Correct_Resampled_Eps.m
%
%---------------------
% May 29, 2015 - A. Pickering
% 11/23/15 - AP - updated
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
% plot eps along sampling path for each realization

eps_prof=nan*ones(length(REsamp.z),REsamp.Nshift);

for whc=1:REsamp.Nshift
    
    clear idt ec tc zc thdir dir_corr Ltc n2c
    
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
        ec=squeeze(REsamp.eps(:,idt(1:2:end),whc));
        tc=squeeze(REsamp.tsamp(:,idt(1:2:end),whc));
        zc=squeeze(REsamp.zsamp(:,idt(1:2:end),whc));
        Lotc=squeeze(REsamp.Lot(:,idt(1:2:end),whc));        
        Ltc=squeeze(REsamp.Lttot(:,idt(1:2:end),whc));
        n2c=squeeze(REsamp.n2(:,idt(1:2:end),whc));
    elseif thdir*dir_corr==-1 % 1st profile is up
        ec=squeeze(REsamp.eps(:,idt(2:2:end),whc));
        tc=squeeze(REsamp.tsamp(:,idt(2:2:end),whc));
        zc=squeeze(REsamp.zsamp(:,idt(2:2:end),whc));
        Lotc=squeeze(REsamp.Lot(:,idt(2:2:end),whc));
        Ltc=squeeze(REsamp.Lttot(:,idt(2:2:end),whc));
        n2c=squeeze(REsamp.n2(:,idt(2:2:end),whc));
    end
    
    % remove 'bad' points
    clear ww
    %   ww=REsamp.Lttot(:,:,whc).*sqrt(REsamp.n2(:,:,whc));
    ww=Ltc.*sqrt(n2c);
%    ww=Lotc.*sqrt(n2c);
    ww=real(ww);
    clear idb
    idb=find(ww>(thresh*REsamp.w_samp));
    ec(idb)=nan;
    
    % save the mean depth profiles
    eps_prof(:,whc)=nanmean(ec,2);
    
    if strcmp(whdir,'up')
        % REsamp.eps is always increasing depth, so we need to flip it for
        % plotting on upward paths
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