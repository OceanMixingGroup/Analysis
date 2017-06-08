function [Epsout,Lmin,Lot,runlmax,Lttot,ptmp]=compute_overturns_discrete(p,t,s,lat,usetemp,minotsize,sigma,runlmin);

%compute_overturns_APcomment.m
%
% %compute_overturns.m, with added comments by AP (trying to figure out
% exactly what is being done etc.) All AP comments preceded by %~ .8 Jan 2015
%
%
% calculate overturns from a variety of instruments
% p,t,s are vectors of pressure, temperature and salinity typically from downcast only
% usetemp set to 1 to use temperature to compute overturns, otherwise uses density
% minotsize is the minimum overturn size to consider (too small and may be noise)
% sigma is noise level for density.
% runlmin is run length minimum.
%
% same as compute_overturns.m but input is salinity instead of conductivity
% becuase sometimes salinity has been separately despiked already
%
%~ OUTPUTS:
%~ -Lttot is 'thorpe scale' (rms displacement over each overturning region).
%~  This is value that is used to compute eps?
%~ -Lmin is minimum overturn resolvable based on noise in density (sigma)
%~ and N2
%~ -runlmax is max run length
%~ -Lot is all the individual thorpe displacements (as opposed to the rms thorpe scale Lttot)
%~ ** 11 Feb - I think Lot is actually the patch size, not the individual
%displacements ***
%~ - ptmp is potential temp?
%
%~ would be nice to return displacements ('dz' I think in this code) and
%sorted profiles for plotting and diagnostics. Also would be nice to return
%a params structure with all the parameters and info.
%
warning off
%% set defaults
if ~exist('usetemp')|isempty(usetemp); usetemp=0; end
if ~exist('lat')|isempty(lat); lat=30; end
if ~exist('minotsize')|isempty(minotsize); minotsize=2; end
if ~exist('sigma')|isempty(sigma); sigma=5e-4; end
if ~exist('runlmin')|isempty(runlmin); runlmin=0; end

%% make potential density and temp at reference depths
% Use one depth if total depth range <1200 m (e.g. fast ctd), but use
% several depth reference levels for a broader range (e.g. shipboard ctd)

if (max(p)-min(p))>1500
    dref=1000; refd=(min(p)+dref/2):dref:max(p);
else
    refd=(min(p)+max(p))/2; dref=(max(p)-min(p))+1;
end

%%
Epsout=NaN*p(:);
Lmin=NaN*t; Lot=NaN*t; Lttot=Lot;


for iref=1:length(refd) %~ loop through different reference densities
    %s = sw_salt(c(:)*10/sw_c3515,t(:),p(:));
    pden = sw_pden(s(:),t(:),p(:),refd(iref));
    ptmp = sw_ptmp(s(:),t(:),p(:),refd(iref));
    
    if usetemp
        V=ptmp;
    else
        V=pden;
    end
    
    %~ sort density profile
    [xx,isort]=sort(pden);
    % smoothed nsq profile
    %if length(t)>600
    %    a=1; b=ones(200,1)/200;
    %    [n2,q,p_ave] = sw_bfrq(nanfilt(b,a,s),nanfilt(b,a,t),p,lat);
    %elseif length(t)>300
    %    a=1; b=ones(100,1)/100;
    %    [n2,q,p_ave] = sw_bfrq(nanfilt(b,a,s),nanfilt(b,a,t),p,lat);
    %else
    
    %~ compute buoynacy freq with sorted profile
    [n2,q,p_ave] = sw_bfrq(s(isort),t(isort),p,lat);
    %end
    
    %~ find good density values (not NaN)
    ig=find(~isnan(pden));
    
    p0=p(:);  %~ oriignal pressure vector
    %~ non-NaN pressure, ptemp, pden, and V
    pg=p(ig); ptmp=ptmp(ig); pden=pden(ig); V=V(ig);
    
    %~ make pg a column vector
    pg=pg(:);
    
    %~ pg is vector of pressures (w/o NaNs)
    %~ V is vector of potential density or temp
    
    %~ not sure what this section is doing...
    %~ find if density is on average increasing or decreasing with depth?
    %~ is this just to distinguish upcast vs downcast?
    sig = sign(nanmedian(diff(V)));
    
    %~ sig*V=??
    %~ sorting by sign and magnitude of dV/dz?
    [tsort,ind]=sort(sig*V);
    %~ why is tsort redefined right away?
    tsort=sig*V;
    psort = pg(ind);
    dz = pg-psort; %~ thorpe displacement?
    
    csdz = cumsum(-dz);
    thresh = 0.0000001;
    
    %~ find starts of overturn regions?
    start = find(csdz(1:end-1)<thresh & csdz(2:end)>=thresh)+1;
    if dz(1)<0
        start = [1;start];
    end;
    %~ find ends of overturn regions?
    stops = find(csdz(1:end-1)>=thresh & csdz(2:end)<thresh)+1;
    %~ make empty arrays for results
    Otnsq = NaN*dz;
    Lt=NaN*dz;
    Lmin0=NaN*pg;
    Lot0=NaN*pg;
    runlmax0=Lmin0;
    R0tot=Lmin0;
    %~ cycle through each region
    for j = 1:length(start);
        %~ find indices of this overturn region
        ind=clip([(start(j)-1):(stops(j)+1)],1,prod(size(dz)));
        %~ get sorted N2 over this region
        indp=find(p_ave>min(pg(ind))&p_ave<max(pg(ind)));
        %~ average N2 of this region
        n2avg=nanmean(n2(indp));
        warning off
        %~ size of region (patch?)
        delz=abs(max(pg(ind))-min(pg(ind)));
        %~ maximum drho over this region
        drhodz=(max(pden(ind))-min(pden(ind)))/delz;
        % run length
        stopnow=0; runl=1;
        ig=find(diff(sign(dz(ind)))==0);
        if ~isempty(ig)
            runl=runl+1;
            ig2=find(diff(ig)==1);
            while stopnow==0
                if isempty(ig2)
                    stopnow=1;
                else
                    ig2=find(diff(ig2)==1);
                    runl=runl+1;
                end
            end
        end
        runlmax0(ind)=runl;
        %~ compute min overturn resolvable based on noise level of density
        %(sigma) and buoyancy freq
        Lmin0(ind)=2*9.8/n2avg*sigma/1027;
        %~ aren't these same as delz and drhodz above?
        Lot0(ind)=(max(pg(ind))-min(pg(ind)));
        drho=(max(pden(ind))-min(pden(ind)));
        %    if (delz>minotsize)&  (length(ind)>(10*(sigma/drho)^2))  % jody's suggestion
        % additional test from Gargett and Garner 08
        Lpos=length(find((V(ind)-sort(V(ind)))>0)); Lneg=length(find((V(ind)-sort(V(ind)))<0));
        R0=min(Lpos/length(ind),Lneg/length(ind));
        
        %~ check if overturn > minimum overturn size and passes run-length
        % test
        if (delz>minotsize)&(delz>(2*9.8/n2avg*sigma/1027))&((max(pden(ind))-min(pden(ind)))>(2*sigma))...
                &runl>runlmin&(max(abs(V(ind)-sort(V(ind))))>2*sigma)&R0>0.2
            Otnsq(ind) = 9.8./mean(pden(ind)).*drhodz;
            temptemp(j)=(max(pg(ind))-min(pg(ind)));
            %~ Lt is rms displacement (thorpe scale)
            Lt(ind)=sqrt(mean(dz(ind).^2));
            R0tot(ind)=R0;
        else
            Otnsq(ind)=NaN; Lmin0(ind)=NaN; Lot0(ind)=NaN; Lt(ind)=NaN; R0tot(ind)=NaN;
        end
    end;
    
    %~ p0 is full-depth pressure vector
    %~ find section of p0 in this reference region
    iz=find(p0>(refd(iref)-dref/2)&p0<=(refd(iref)+dref/2));
    
    [xxx,iun]=unique(pg); Lt=Lt(:);
    
    %~ interpolate Lt back to full depth vector
    %~ Lt and Lttot are rms thrope scales over regions
    Lttot(iz)=interp1(pg(iun),Lt(iun),p0(iz));
    %~ compute eps and inteprolate to full depth vector
    Epsout(iz) = interp1(pg(iun),0.64*Lt(iun).^2.*sqrt(Otnsq(iun)).^3,p0(iz));
    %~ interp Lmin, Lot, runlmax to full depth vector
    Lmin(iz)=interp1(pg(iun),Lmin0(iun),p0(iz));
    %~ Lot is patch size?
    Lot(iz)=interp1(pg(iun),Lot0(iun),p0(iz));
    runlmax(iz)=interp1(pg(iun),runlmax0(iun),p0(iz));
    
    %~ fill in NaN eps values
    Epsout(isnan(Epsout))=1e-11;
    
    
    
end

