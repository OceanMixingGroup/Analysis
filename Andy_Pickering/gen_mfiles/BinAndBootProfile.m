function [bb, zbins, hm]=BinAndBootProfile(xin,zin,interval,nboot,plotit)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% [xout bb hm]=BinAndBootProfiles(xin,zin,interval,nboot,plotit)
%
% Average data from profile of some variable 'xin' in depth bins of size 'interval',
% and then compute mean and 95% bootstrap confidence limits. Option to plot.
%
% Can be also be more than 1 profile concatenated together.
%
%
% INPUT:
% xin      - vector of data to bootstrap [M X 1]
% zin      - depths corresponding to xin [M X 1]
% interval - interval to bootstrap over [1X1]
% nboot    - # times to bootstrap [1 X 1]
% plotit   - option to plot (0,1)
%
% OUTPUT:
% xout - the bootsrapped data - length=ceil(nanmax(zin)/interval)
% bb   - [low_conf mean upper_conf]
% hm   - plot handle
%
% Depends:
% - boot_v5.m
%
%----------------------------
% A. Pickering - July 20 2015 - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

addpath /Users/Andy/Cruises_Research/mixingsoftware/general/

nbins = ceil(nanmax(zin)/interval);
bins  = 0:interval:(nbins*interval);
NN    = nan(ones(1,nbins));

% mid-points of bins
zbins = interval/2 : interval : (nbins*(interval));

% Make empty arrays
xout     = nan*ones(length(bins)-1,1);
boot_low = nan*ones(length(bins)-1,1);
boot_up  = nan*ones(length(bins)-1,1);

% loop over depth bins
for whbin=1:nbins-1
    
    clear idz b dataToboot
    
    % find data in this depth bin
    idz = find(zin>bins(whbin) & zin<bins(whbin+1));
    %xout(idz,:) = repmat(nanmean(xin(idz,:)),length(idz),1);
    
    % bootstrap data for this bin
    dataToboot = xin(idz);
    % get rid of Nans
    dataToboot = dataToboot(~isnan(dataToboot));
    [b] = boot_v5(dataToboot,nboot);
    
    % bootstrap mean
    xout(whbin) = b(2);
    % lower conf. limit
    boot_low(whbin) = b(1);
    % upper conf. limit
    boot_up(whbin) = b(3);
    
end

bb = [boot_low xout boot_up];

% Plotting w/ shading of confidence limits
if plotit==1
    
    hm = plot(log10(bb(:,2)),zbins,'kd','linewidth',4) ;
    hold on
    
    for iz = 1:length(zbins)
        line([log10(bb(iz,[1 3]))],[zbins(iz) zbins(iz)],'color','k')
    end
    
    axis ij
    grid on
else
    hm=[];
end % plotit
%%