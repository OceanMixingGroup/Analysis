%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% SampleSynthFieldsOT.m
%
% Make a synthetic temperature field with stable (no overturns) sinusoidal
% oscillations, and resample with simulated MP or CTD at different speeds
% to see if false overturns are introduced.
%
%
% 9 Feb 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3/Tchain3_RecomputedEps_MinOT_50.mat')

% make a synthetic temperature field

t=0 : 2/60/24 : 1.5 ;

% ~~ Use a linear temperature profile
z=0:10:2000;
tmean=linspace(2,8,length(z));
tmean=flipud(tmean');

% ~~ Use actual mean temp profile from T-chain
% z=xx2.z;
% tmean=nanmean(xx2.T,2);
% z=z(~isnan(tmean));
% tmean=tmean(~isnan(tmean));

[T,Z]=meshgrid(t,z);

% ~~ add a sinusoidal variability
Tsyn=12 % period of oscillations (hr)
om=2*pi*24/Tsyn; % frequency of oscillations
A=0.2
A=0.7
temp=nan*ones(size(T));
temp=repmat(tmean,1,length(t))+ A.*sin(om.*T);

figure(1);clf
agutwocolumn(0.6)
wysiwyg
ezpc(T,Z,temp)
hold on
%contour(T,Z,temp,[10:30],'k')
contour(T,Z,temp,tmean(1:10:end),'k')
cb=colorbar
cb.Label.String='temp'
xlabel('Time [days]')
ylabel('Depth [m] ')


%% Now simulate sampling

%~ Sampling parameters
clear SP
SP.z_range=[z(1) z(end)];
SP.dz_samp=10;
SP.w_samp=0.045;
SP.t_start=0.0;     % start time
SP.tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
SP.time_range=1.5;

%~ the 'true' data to sample
data_real=temp;
t_real=t;
z_real=z;

addpath /Users/Andy/Cruises_Research/SimProfiler/

REsamp=ResampleFieldGeneral(data_real,t_real,z_real,SP)

%% add overturns to resampled data

Params.lat=20.5;
Params.usetemp=0;
Params.minotsize=10;
Params.sigma=1e-5;
Params.runlmin=0;
Params.plotit=0;

addpath /Users/Andy/Cruises_Research/SimProfiler/

REsamp=AddOverturnsToREsamp(REsamp,Params)

%%

makemovie=0
if makemovie==1
    clear nFrames mov whframe
    nFrames=length(1:5:REsamp.Nshift);
    % Preallocate movie structure.
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    whframe=1
end
%
tm=nanmean(REsamp.data_real,2);
dd=12

for whc=1:5:REsamp.Nshift
    figure(1);clf
    agutwocolumn(0.7)
    wysiwyg
    ezpc(REsamp.treal,REsamp.zreal,REsamp.data_real);
    hold on
    contour(REsamp.treal,REsamp.zreal,REsamp.data_real,tm(1:dd:end),'k');
    plot(REsamp.tsamp(:,:,whc),REsamp.zsamp(:,:,whc),'w--');
    %  plot(REsamp.tsamp(:,:,1),REsamp.zsamp(:,:,1),'w')
    contour(REsamp.tgrid(whc,:),REsamp.z,REsamp.data_resamp(:,:,whc),tm(1:dd:end),'w');
    vline(REsamp.tgrid(whc,:),'--');
    ylabel('Depth [m]')
    xlabel('Time [days]')
    cb=colorbar
    
    Nprof=size(REsamp.data_resamp,2);
    zmat=REsamp.zsamp(:,:,whc);
    tmat=REsamp.tsamp(:,:,whc);
    for whp=1:Nprof
        clear t2 dtdz idot
        if nanmean(diff(zmat(:,whp)))<0
            t2=flipud(tmat(:,whp));
        else
            t2=tmat(:,whp);
        end
        dtdz=diffs(REsamp.data_resamp(:,whp,whc));
        idot=find(dtdz>0);
        plot(tmat(:,whp),zmat(:,whp),'y');
        plot(t2(idot),z(idot),'ys');
        shg
    end
    pause(0.15)
    
    if makemovie==1
        mov(whframe)=getframe(gcf);
        whframe=whframe+1;
    end
end

if makemovie==1
    savdir='/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/NotesOverturnBiases/'
    fname=fullfile(savdir,['SyntheticExample.avi']);
    movie2avi(mov, fname, 'compression', 'None','fps',1);
end

%%