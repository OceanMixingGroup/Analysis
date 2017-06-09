% Script to calibrate sensors:
global head data cal
modify_header_eq08
%%
% calibrate is called as 
% calibrate('series','method',{'filter1','filter2',...})
% where series is the series to calibrate ( 'T1', 'TP2', 'S1', 'UC', etc.)
% and method is the method to use ('T','TP','S','UC', etc.)
% filter is 'Axxx' where A is h,l,n for highpass, lowpass or notch filter
% and xxx is the cutoff frequency.  If A='n', then Axxx='n20-25' would notch
% out the frequencies between 20 and 25 Hz.
% NOTE: series could be 'temp' if data.TEMP coef.TEMP and irep.TEMP all exist
% prior to calling calibrate.

% bad=0;
% if q.script.num~=2145 && q.script.num~=2146
%     [data,bad]=issync(data,head);
% end
if strfind(head.module_num(head.sensor_index.TP,:),'04-01a') 
    head.coef.TP=[1 0.09495 0 0 1];
elseif strfind(head.module_num(head.sensor_index.TP,:),'04-02a') 
    head.coef.TP=[1 0.086257 0 0 1];
elseif strfind(head.module_num(head.sensor_index.TP,:),'04-03a') 
%     head.coef.TP=[1 1 -0.49696 -1.5154 1]; % correct for 04-03a 
% however, 04-05a circuit is installed in CHAM04-03,
% while header reads 04-03a
% so these coefficiants are for 04-05a:
    head.coef.TP=[1 0.094198 0 0 1];
elseif strfind(head.module_num(head.sensor_index.TP,:),'04-04a') 
    head.coef.TP=[1 0.083988 0 0 1];
end
switch q.script.num
    case {1861,1862,1871,1872,2145,2146,2150,2191,2617,2618,2619,2620}
        bad=0;
    case 2473
        bad=0;
        fields=fieldnames(data);
        for jj=1:length(fields)
            data.(char(fields(jj)))=data.(char(fields(jj)))...
                (515*head.irep.(char(fields(jj)))-head.irep.(char(fields(jj)))+1:end);
        end
    otherwise
        [data,bad]=issync(data,head);
end
        
cg=[1861,1862,1871,1872,2145,2146,2150,2191,2355:2358,2617:2620];
for ii=cg
    if q.script.num==ii
        data.S1(data.S1<-3)=NaN;
        data.S2(data.S2<-3)=NaN;
        data.TP(data.S2<-1)=NaN;
        data.SCAT(data.SCAT<0)=NaN;
        for kk=1:3
            data.MHC=deglitch(data.MHC,200,2);
            data.T=deglitch(data.T,200,2);
            data.P=deglitch(data.P,200,2);
            if isfield(data,'MHT')
                data.MHT=deglitch(data.MHT,200,2);
            end
        end
        data.P=fillgap(data.P);
        data.T=fillgap(data.T);
        data.MHC=fillgap(data.MHC);
        data.S1=fillgap(data.S1);
        data.S2=fillgap(data.S2);
        data.TP=fillgap(data.TP);
        data.SCAT=fillgap(data.SCAT);
        if isfield(data,'MHT')
            data.MHT=fillgap(data.MHT);
        end
    end
end

% for kk=1:2
%     data.MHC=deglitch(data.MHC,round(length(data.MHC)/7),2);
%     data.T=deglitch(data.T,round(length(data.MHC)/7),2);
%     data.P=deglitch(data.P,round(length(data.MHC)/7),2);
%     if isfield(data,'MHT')
%         data.MHT=deglitch(data.MHT,round(length(data.MHC)/7),3);
%     end
% end
% data.P=fillgap(data.P);
% data.T=fillgap(data.T);
% data.MHC=fillgap(data.MHC);
% if isfield(data,'MHT')
%     data.MHT=fillgap(data.MHT);
% end

% if q.script.num==1861
%     calibrate('p','p','l0.1');
% end
if bad==1
   disp('no sync');
    return;
end
if length(data.P)<200
   disp('no sync or too short file');
    bad=1;
    return;
end
calibrate('p','p','l2');
calibrate('az','az')
calibrate('ax','ax')
calibrate('ay','ay')
% if q.script.num==1861
%     calibrate('p','p','l0.1');
% end
if (abs(max(cal.P)-min(cal.P))<4 || any(cal.P<-3));
   bad=1;
   disp('profile < 4m or negative P');
   disp(['min pressure = ' num2str(min(cal.P))]);
   return;
end;
% determine if it's down or up cast and set the flag
if nanmean(diff(cal.P))<0
   head.direction='u';
else
   head.direction='d';
end
% run some script or function which selects the appropriate depth range and
% places the indices into q.mini and q.maxi
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% determine_depth_range is one possibility...
if head.direction=='d'
  [q.mini,q.maxi,head.got_bottom]=determine_depth_range2(0.5);

  % now select only the data within that depth range
  len=select_depth_range(q.mini,q.maxi);
  [data.P,mini,maxi]=extrapolate_depth_range(data.P,min(5000,(length(cal.P)-2)));
  % extrapolate_depth_range flips the ends of p over itself before calibrating so
  % that starting and ending transients are eliminated.
else
   mini=1; maxi=length(data.P);
   head.got_bottom='n';
end
if maxi-mini<200
    disp('profile < 4m');
    bad=1;
    return;
end
calibrate('p','p','l2') 
calibrate('p','fallspd','l.5') 
data.P=data.P(mini:maxi);
cal.P=cal.P(mini:maxi);
if (abs(max(cal.P)-min(cal.P))<3);
   bad=1;
   disp('profile < 4m');
   return;
end;
calibrate('ay','tilt','l.8')
calibrate('ax','tilt','l.8')
temp=filter_series(cal.AZ,102.4*head.irep.AZ,'h2'); % this was changed from .2 to 2
cal.AZ2=temp.*temp;
head.irep.AZ2=head.irep.AZ;
cal.FALLSPD=cal.FALLSPD(mini:maxi);
q.fspd=mean(cal.FALLSPD);
% calibrate('s1','s',{'h.4'})
% calibrate('s2','s',{'h.4'})
calibrate('s1','s',{'h0.4','l50'})
calibrate('s2','s',{'h0.4','l50'})
% calibrate('s1','s',{'h0.4','l10'})
% calibrate('s2','s',{'h0.4','l10'})
% Determine the cutoff frequency for w:
% this should be a cutoff frequency at 3m.
freq=num2str(2*q.fspd/100/3) ;
%calibrate('w','w',{['h' freq]})
%calibrate('w','volts','h0.01')
spd=q.fspd;
%sp=1/head.coef.W(2);
rho=1.024;
%cal1.W=cal.W./(2*rho*spd.*sp);
%calibrate('w','volts',{['h' freq]})
%cal2.W=cal.W./(2*rho*spd.*sp);
calibrate('t','t')
cal.T(cal.T>26)=NaN;cal.T(cal.T<11)=NaN;
if isfield(data,'MHT')
    calibrate('mht','t')
else
    cal.MHT=cal.T;
    head.irep.MHT=head.irep.T;
    cal.MHT(cal.MHT>26)=NaN;cal.MHT(cal.MHT<11)=NaN;
end

% calibrate('tp','volts',{'h1','l10'})
% temp=cal.TP;
calibrate('tp','tp')
% cal.TP=calibrate_tp(data.TP,head.coef.TP,data.T,head.coef.T,cal.FALLSPD);
% data.TP=temp;

calibrate('MHC','c','l20')
cal.C=cal.MHC;head.irep.C=head.irep.MHC;
cal.C(cal.C<2)=NaN;cal.C(cal.C>7)=NaN;
if q.script.num==544
    cal.C(cal.P<92.1)=NaN;
end
calibrate('scat','volts')

cond=cal.C(1:head.irep.C:length(cal.C));
press=cal.P;
% if isfield(data,'MHT')
%     temp=cal.MHT(1:head.irep.MHT:end);
%     cal.S=sw_salt(10*cond/sw_c3515,temp,press); % convert to mmho/cm 1st
%     head.irep.S=head.irep.P;
%     calc_theta('theta','s','mht','p');
%     calc_sigma('sigma','s','mht','p');
% else
temp=cal.T(1:head.irep.T:end);
cal.S=sw_salt(10*cond/sw_c3515,temp,press); % convert to mmho/cm 1st
% cal.S=deglitch(cal.S,100,4);
% cal.S=deglitch(cal.S,100,3);
% cal.S=deglitch(cal.S,100,3);
% cal.S=fillgap(cal.S);
head.irep.S=head.irep.P;
calc_theta('theta','s','t','p');
calc_sigma('sigma','s','t','p');
cal.SIGMA=fillgap(cal.SIGMA);
% end
%cal.SIGTH=sw_pden(salinity,cal.T,cal.P,0)-1000head;
%head.irep.SIGTH=head.irep.P;

% calculate SIGMA_ORDER
inds=calc_order('sigma','P');
% calculate THETA_ORDER
inds=calc_order('theta','P');
% calculate the mean temperature gradient on small scales:
cal.DTDZ=diff(cal.THETA_ORDER)./diff(cal.P);
cal.DTDZ(length(cal.DTDZ)+1)=cal.DTDZ(length(cal.DTDZ));
head.irep.DTDZ=1;
% calculate the mean density gradient on small scales:
cal.DRHODZ=diff(cal.SIGMA_ORDER)./diff(cal.P);
cal.DRHODZ(length(cal.DRHODZ)+1)=cal.DRHODZ(length(cal.DRHODZ));
head.irep.DRHODZ=1;
g=9.81;
cal.SIGMA=cal.SIGMA-1000;
cal.SIGMA_ORDER=cal.SIGMA_ORDER-1000;
rhoav=mean(cal.SIGMA(1:1:length(cal.SIGMA)-1))+1000;
cal.N2=(g/rhoav).*diff(cal.SIGMA_ORDER)./diff(cal.P);
cal.N2(length(cal.N2)+1)=cal.N2(length(cal.N2));
head.irep.N2=head.irep.P;
cal.VARLT=(cal.THORPE_SIGMA-mean(cal.THORPE_SIGMA)).^2; %variance Thorpe scale 
cal.VARLT_THETA=(cal.THORPE_THETA-mean(cal.THORPE_THETA)).^2; %variance Thorpe scale for temperature 
head.irep.VARLT=head.irep.THORPE_SIGMA;
head.irep.VARLT_THETA=head.irep.THORPE_SIGMA;
cal.VARAZ=(cal.AZ-mean(cal.AZ)).^2; %variance of AZ
head.irep.VARAZ=head.irep.AZ;

if ~exist('cast','var')
    cast=q.script.num;
end
tag_file_eq08;
for jjjj=1:size(tag.c,1)
    if cast==floor(tag.c(jjjj,1)) && cast==floor(tag.c(jjjj,2))
        bad_c=find(cal.P>=round((tag.c(jjjj,1)-floor(tag.c(jjjj,1)))*1000)-0.6 & cal.P<=round((tag.c(jjjj,2)-floor(tag.c(jjjj,2)))*1000)+0.6);
    elseif cast==floor(tag.c(jjjj,1)) && cast<floor(tag.c(jjjj,2))
        bad_c=find(cal.P>=round((tag.c(jjjj,1)-floor(tag.c(jjjj,1)))*1000)-0.6);
    elseif cast>floor(tag.c(jjjj,1)) && cast<floor(tag.c(jjjj,2))
        bad_c=find(cal.P);
    elseif cast>floor(tag.c(jjjj,1)) && cast==floor(tag.c(jjjj,2))
        bad_c=find(cal.P<=round((tag.c(jjjj,2)-floor(tag.c(jjjj,2)))*1000)+0.6);
    else
        bad_c=[];
    end
    if(~isempty(bad_c))
        bd={'C','S','SIGMA','SIGMA_ORDER','N2','DRHODZ','VARLT','VARLT_THETA'};
        for iiii=1:length(bd)
            eval(['irep=head.irep.' char(bd(iiii)) ';']);
            eval(['cal.' char(bd(iiii)) '(bad_c(1)*irep:bad_c(end)*irep)=NaN;']);
        end
    end
end
for jjjj=1:size(tag.t,1)
    if cast==floor(tag.t(jjjj,1)) && cast==floor(tag.t(jjjj,2))
        bad_t=find(cal.P>=round((tag.t(jjjj,1)-floor(tag.t(jjjj,1)))*1000)-0.6 & cal.P<=round((tag.t(jjjj,2)-floor(tag.t(jjjj,2)))*1000)+0.6);
    elseif cast==floor(tag.t(jjjj,1)) && cast<floor(tag.t(jjjj,2))
        bad_t=find(cal.P>=round((tag.t(jjjj,1)-floor(tag.t(jjjj,1)))*1000)-0.6);
    elseif cast>floor(tag.t(jjjj,1)) && cast<floor(tag.t(jjjj,2))
        bad_t=find(cal.P);
    elseif cast>floor(tag.t(jjjj,1)) && cast==floor(tag.t(jjjj,2))
        bad_t=find(cal.P<=round((tag.t(jjjj,2)-floor(tag.t(jjjj,2)))*1000)+0.6);
    else
        bad_t=[];
    end
    if(~isempty(bad_t))
        bd={'T','TP','S','THETA','SIGMA','SIGMA_ORDER','N2','DTDZ','DRHODZ','VARLT','VARLT_THETA'};
        for iiii=1:length(bd)
            eval(['irep=head.irep.' char(bd(iiii)) ';']);
            eval(['cal.' char(bd(iiii)) '(bad_t(1)*irep:bad_t(end)*irep)=NaN;']);
        end
    end
end
for jjjj=1:size(tag.mht,1)
    if cast==floor(tag.mht(jjjj,1)) && cast==floor(tag.mht(jjjj,2))
        bad_t=find(cal.P>=round((tag.mht(jjjj,1)-floor(tag.mht(jjjj,1)))*1000)-0.6 & cal.P<=round((tag.mht(jjjj,2)-floor(tag.mht(jjjj,2)))*1000)+0.6);
    elseif cast==floor(tag.mht(jjjj,1)) && cast<floor(tag.mht(jjjj,2))
        bad_t=find(cal.P>=round((tag.mht(jjjj,1)-floor(tag.mht(jjjj,1)))*1000)-0.6);
    elseif cast>floor(tag.mht(jjjj,1)) && cast<floor(tag.mht(jjjj,2))
        bad_t=find(cal.P);
    elseif cast>floor(tag.mht(jjjj,1)) && cast==floor(tag.mht(jjjj,2))
        bad_t=find(cal.P<=round((tag.mht(jjjj,2)-floor(tag.mht(jjjj,2)))*1000)+0.6);
    else
        bad_t=[];
    end
    if(~isempty(bad_t))
        bd={'MHT'};
        for iiii=1:length(bd)
            eval(['irep=head.irep.' char(bd(iiii)) ';']);
            eval(['cal.' char(bd(iiii)) '(bad_t(1)*irep:bad_t(end)*irep)=NaN;']);
        end
    end
end
for jjjj=1:size(tag.sig,1)
    if cast==floor(tag.sig(jjjj,1)) && cast==floor(tag.sig(jjjj,2))
        bad_t=find(cal.P>=round((tag.sig(jjjj,1)-floor(tag.sig(jjjj,1)))*1000)-0.6 & cal.P<=round((tag.sig(jjjj,2)-floor(tag.sig(jjjj,2)))*1000)+0.6);
    elseif cast==floor(tag.sig(jjjj,1)) && cast<floor(tag.sig(jjjj,2))
        bad_t=find(cal.P>=round((tag.sig(jjjj,1)-floor(tag.sig(jjjj,1)))*1000)-0.6);
    elseif cast>floor(tag.sig(jjjj,1)) && cast<floor(tag.sig(jjjj,2))
        bad_t=find(cal.P);
    elseif cast>floor(tag.sig(jjjj,1)) && cast==floor(tag.sig(jjjj,2))
        bad_t=find(cal.P<=round((tag.sig(jjjj,2)-floor(tag.sig(jjjj,2)))*1000)+0.6);
    else
        bad_t=[];
    end
    if(~isempty(bad_t))
        bd={'THETA','SIGMA','N2','DTDZ','DRHODZ','VARLT','VARLT_THETA'};
        for iiii=1:length(bd)
            eval(['irep=head.irep.' char(bd(iiii)) ';']);
            eval(['cal.' char(bd(iiii)) '(bad_t(1)*irep:bad_t(end)*irep)=NaN;']);
        end
    end
end
for jjjj=1:size(tag.n2,1)
    if cast==floor(tag.n2(jjjj,1)) && cast==floor(tag.n2(jjjj,2))
        bad_t=find(cal.P>=round((tag.n2(jjjj,1)-floor(tag.n2(jjjj,1)))*1000)-0.6 & cal.P<=round((tag.n2(jjjj,2)-floor(tag.n2(jjjj,2)))*1000)+0.6);
    elseif cast==floor(tag.n2(jjjj,1)) && cast<floor(tag.n2(jjjj,2))
        bad_t=find(cal.P>=round((tag.n2(jjjj,1)-floor(tag.n2(jjjj,1)))*1000)-0.6);
    elseif cast>floor(tag.n2(jjjj,1)) && cast<floor(tag.n2(jjjj,2))
        bad_t=find(cal.P);
    elseif cast>floor(tag.n2(jjjj,1)) && cast==floor(tag.n2(jjjj,2))
        bad_t=find(cal.P<=round((tag.n2(jjjj,2)-floor(tag.n2(jjjj,2)))*1000)+0.6);
    else
        bad_t=[];
    end
    if(~isempty(bad_t))
        bd={'N2'};
        for iiii=1:length(bd)
            eval(['irep=head.irep.' char(bd(iiii)) ';']);
            eval(['cal.' char(bd(iiii)) '(bad_t(1)*irep:bad_t(end)*irep)=NaN;']);
        end
    end
end
