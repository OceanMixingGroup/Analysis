% cali_eq14.m
%
% adapted from cali_dn11a by Sally Warner, March 2015
%
% many changes made to variable names  when Pavan took over the chameleon
%   -   conductivity: C was changed to COND
%   -   temperature: T and T2 were changed to T1 and T2
%   -   salinity: S changed to SAL
%
% Comments from old cali_realtime:
%
% calibrate is called as 
% calibrate('series','method',{'filter1','filter2',...})
% where series is the series to calibrate ( 'T1', 'TP2', 'S1', 'UC', etc.)
% and method is the method to use ('T','TP','S','UC', etc.)
% filter is 'Axxx' where A is h,l,n for highpass, lowpass or notch filter
% and xxx is the cutoff frequency.  If A='n', then Axxx='n20-25' would notch
% out the frequencies between 20 and 25 Hz.
% NOTE: series could be 'temp' if data.TEMP coef.TEMP and irep.TEMP all exist
% prior to calling calibrate.

%% Script to calibrate sensors:
global head data cal

%% modify the header coefficients for conductivity (and other variablies)

modify_header_eq14

%% determine if the cast is good

cast=q.script.num;

bad=0;

% make sure the chameleon has been synced with the computer
[data,bad]=issync(data,head);
if bad==1
    disp('no sync');
    return;
end

% make sure length of file is long enough
if length(data.P)<200
    disp('no sync or too short file');
    bad=1;
    return;
end


%% calibrate  pressure and accelerations

% calibrate pressure and vertical accelerations and check to make sure
% profile is long enough and that pressure is right sign
% also calibrate horizontal accelerations
calibrate('p','p','l2');
calibrate('az','az')
calibrate('ax','ax')
calibrate('ay','ay')
if (abs(max(cal.P)-min(cal.P))<4 || any(cal.P<-3))
   bad=1;
   disp('profile < 4m or negative P');
   disp(['min pressure = ' num2str(min(cal.P))]);
   return
end;

% determine if it's down or up cast and set the flag (sasha's comment)
if nanmean(diff(cal.P))<0
   head.direction='u';
else
   head.direction='d';
end
% run some script or function which selects the appropriate depth range and
% places the indices into q.mini and q.maxi (sasha's comment)

% determine_depth_range is one possibility... (sasha's comment)
if head.direction=='d'
  [q.mini,q.maxi,head.got_bottom]=determine_depth_range2(0.5);
  % note in the readtime cali file it was determine_depth_range2(3)

  % now select only the data within that depth range (sasha's comment)
  len=select_depth_range(q.mini,q.maxi);
  [data.P,mini,maxi]=extrapolate_depth_range(data.P,min(5000,(length(cal.P)-2)));
  % extrapolate_depth_range flips the ends of p over itself before calibrating so
  % that starting and ending transients are eliminated. (sasha's comment)
else
   mini=1; maxi=length(data.P);
   head.got_bottom='n';
end

% check to make sure the cast is long enough
if maxi-mini<200
    disp('profile < 4m');
    bad=1;
    return;
end

%%%% PRESSURE %%%%
calibrate('p','p','l2') 
calibrate('p','fallspd','l.5') 
data.P=data.P(mini:maxi);
cal.P=cal.P(mini:maxi);
if (abs(max(cal.P)-min(cal.P))<3);
   bad=1;
   disp('profile < 4m');
   return;
end;

%%%% ACCELERATIONS/TILTS %%%%
calibrate('ay','tilt','l.8')
calibrate('ax','tilt','l.8')
temp=filter_series(cal.AZ,102.4*head.irep.AZ,'h2'); % this was changed from .2 to 2 (sasha's comment)
cal.AZ2=temp.*temp;
head.irep.AZ2=head.irep.AZ;

%%%% FALL SPEED %%%%
cal.FALLSPD=cal.FALLSPD(mini:maxi);
q.fspd=mean(cal.FALLSPD);

% Determine the cutoff frequency for w: (sasha's comment)
% this should be a cutoff frequency at 3m. (sasha's comment)
% (sjw) the following, except freq, spd and rho were all naned out. And it
% doesn't look like any of them are used...
freq=num2str(2*q.fspd/100/3) ;
%calibrate('w','w',{['h' freq]})
%calibrate('w','volts','h0.01')
spd=q.fspd;
%sp=1/head.coef.W(2);
rho=1.024;
%cal1.W=cal.W./(2*rho*spd.*sp);
%calibrate('w','volts',{['h' freq]})
%cal2.W=cal.W./(2*rho*spd.*sp);


%% calibrate other varibles 

%%%% SHEAR %%%%
calibrate('s1','s',{'h0.4'})
calibrate('s2','s',{'h0.4'})
% was {'h.4','l20'} in realtime version and older versions


%%%% TEMPERATURE %%%%
calibrate('t1','t')  %sjw: changed "t" to "t1"
calibrate('t2','t','l20')


%%%% TEMPERATURE DIFFERENTIAL %%%%
if head.coef.TP(2)==0
    head.coef.TP(2)=0.1;
end
% calibrate('tp','tp') 
% Sasha had been using calibrate('tp','tp') to calibrate TP,
% but I get an error message in calibrate and even if I tweak the calibrate
% function, the result doesn't make sense: I think it uses the temperature 
% (T) coefficients instead of TP. Went back to using the three lines of code that
% were in cali_realtime_v3 which gives a desirable looking result (sjw 3/15)
calibrate('tp','volts',{'h1','l15'});
data.TP=cal.TP;
cal.TP=calibrate_tp(data.TP,head.coef.TP,data.T1,head.coef.T1,cal.FALLSPD);


%%%% CONDUCTIVITY %%%
% MHC has been changed to COND
calibrate('cond','c','l20')
cal.COND(cal.COND<=0)=NaN;

% sasha wrote a long file that corrects the coefficients for the
% conductivity data AGAIN. I don't understand why there are a whole new set
% of coefficients when I thought they were all fixed in modify_header.
% I have not written a new version of fix_MHC, nor do I plan to run this
% function: (commented out next line)
% fix_MHC_dn11a


%%%% SCATTEROMETER %%%
calibrate('scat','volts')


%%%% SALINITY %%%
% S changed to SAL by Pavan
% first need conductiity, pressure and temp
cond=cal.COND(1:head.irep.COND:length(cal.COND));
press=cal.P;
temp=cal.T1(1:head.irep.T1:end);
% use T2 instead of T1 if needed
% if (cast>=1026 && cast<=1217) || (cast>=1695 && cast<=1829) 
%     temp=cal.T2(1:head.irep.T2:end);
% end
% use the seawater toolbox to calculate salinity
cal.SAL=sw_salt(10*cond/sw_c3515,temp,press); % convert to mmho/cm 1st
head.irep.SAL=head.irep.P;


%%%% THETA AND SIGMA %%%
calc_theta('theta','sal','t1','p');
calc_sigma('sigma','sal','t1','p');

% use T2 instead of T1 if needed
% make sure these numbers are set to nan out only T1 but not T2 and THETA
% in the tag file
if cast==2282 ||  cast==2391  || cast==2762 || cast==2953
    calc_theta('theta','sal','t2','p');
    calc_sigma('sigma','sal','t2','p');
end

cal.SIGMA=fillgap(cal.SIGMA);


%% REORDER THE DATA FOR CALCULATING THORPE SCALES ETC
% (sjw 2015) I understand why the data should be reordered for calculating
% dT/dz, Nsq and thorpe scales. However, 

% calculate SIGMA_ORDER
% note: the function calc_order('sigma','P') calculates: THORPE_SIGMA and
% SIGMA_ORDER, and gives the ordered index for density (aka sigma, aka rho)
inds=calc_order('sigma','P');

% then reorder theta with respect to the reordered density index
cal.THETA_RHOORDER=cal.THETA(inds);
head.irep.THETA_RHOORDER=head.irep.THETA;

% then reorder salinity with respect to the reordered density index
cal.SAL_RHOORDER=cal.SAL(inds);
head.irep.SAL_RHOORDER=head.irep.SAL;

% (sjw 3/2015) the following bit was already NaNed out. These seem like
% useful variables... should they be unNaNed or are they calculated
% elsewhere? YES. They are calculated elsewhere in run_eq14 dT/dz, drho/dz
% and N2 are all calculated from 1m averaged data. DO NOT need to be
% calibrated here.
%%%%
% % calculate THETA_ORDER
% inds=calc_order('theta','P');
% % calculate the mean temperature gradient on small scales:
% cal.DTDZ=diff(cal.THETA_ORDER)./diff(cal.P);
% cal.DTDZ(length(cal.DTDZ)+1)=cal.DTDZ(length(cal.DTDZ));
% head.irep.DTDZ=1;
% % calculate the mean density gradient on small scales:
% cal.DRHODZ=diff(cal.SIGMA_ORDER)./diff(cal.P);
% cal.DRHODZ(length(cal.DRHODZ)+1)=cal.DRHODZ(length(cal.DRHODZ));
% head.irep.DRHODZ=1;
% g=9.81;
%%%%

% subtract 1000 from sigma and sigma_order
cal.SIGMA=cal.SIGMA-1000;
cal.SIGMA_ORDER=cal.SIGMA_ORDER-1000;

% (sjw 3/2015) The following calculation of density and N2 was NaNed out
% already
%%%%
% rhoav=mean(cal.SIGMA(1:1:length(cal.SIGMA)-1))+1000;
% cal.N2=(g/rhoav).*diff(cal.SIGMA_ORDER)./diff(cal.P);
% cal.N2(length(cal.N2)+1)=cal.N2(length(cal.N2));
% head.irep.N2=head.irep.P;
%%%%

% calculate variance Thorpe scales and variance of AZ
cal.VARLT=(cal.THORPE_SIGMA-mean(cal.THORPE_SIGMA)).^2; %variance Thorpe scale 
head.irep.VARLT=head.irep.THORPE_SIGMA;
% cal.VARLT_THETA=(cal.THORPE_THETA-mean(cal.THORPE_THETA)).^2; %variance Thorpe scale for temperature 
% head.irep.VARLT_THETA=head.irep.THORPE_SIGMA;
cal.VARAZ=(cal.AZ-mean(cal.AZ)).^2; %variance of AZ
head.irep.VARAZ=head.irep.AZ;

