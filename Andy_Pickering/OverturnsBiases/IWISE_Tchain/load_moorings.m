% script to load in all mooring data.

try
    delete(h)
end
% the following is needed for interp1
warning('off','MATLAB:interp1:NaNinY')

convert_from_cnv=0; % either convert from cnv or read in the matfiles

load mooring_config
config.root_dir='../data/';
config.processed_dir='../data/processed/';
config.raw_path='raw/';
config.delt=2/60/24;
config.mooring_paths={'T1','T2','T3','T4','A1'};
config.sensor_types={'sbe56','sbe37','rbr','chipod_no_data_yet'};
    h=waitbar(0,'loading sensors')
for a=3 % %:4 %1:5 % loop through each mooring 1 to 4 (or 5)
    if a<5
        config.time=datenum('june-11-2011 12:00:00'):config.delt:datenum('aug-25-2011');
        config.time=datenum('june-11-2011 12:00:00'):config.delt:datenum('sep-05-2011');
    else
        config.time=datenum('june-11-2011 12:00:00'):config.delt:datenum('aug-05-2011 00:18:00');
    end
    len_t=length(config.time);
    n_sensors=length(mooring{a}.config.sensornums);  
    init=NaN*ones(n_sensors,len_t);
    mooring{a}.time=config.time;
    mooring{a}.T=init;
    mooring{a}.Tvar=init;
    mooring{a}.P=init;
    for b=1:(n_sensors)
        waitbar(b/n_sensors,h)
        % load data from raw
        if convert_from_cnv
            serial_num=mooring{a}.config.sensornums(b);
            clear dat
        switch mooring{a}.config.sensor_type(b)
            case 1 % SBE-56
                fname=dir([config.root_dir config.raw_path config.mooring_paths{a} '/sbe56/SBE056*' num2str(serial_num) '_*.cnv'])
                dat=load_sb56_cnv2mat([config.root_dir config.raw_path config.mooring_paths{a} '/sbe56/'],fname.name);
            case 2 % SBE-37
                fname=dir([config.root_dir config.raw_path config.mooring_paths{a} '/sbe37/SBE37*' num2str(serial_num) '.asc'])
                if length(fname)==0
                    fname=dir([config.root_dir config.raw_path config.mooring_paths{a} '/sbe37/SBE37*' num2str(serial_num) '*.cnv'])
                end
                dat=read_sbe_jn([config.root_dir config.raw_path config.mooring_paths{a} '/sbe37/' fname.name]);    
            case 3 % RBR
                fname=dir([config.root_dir config.raw_path config.mooring_paths{a} '/rbrs/*' num2str(serial_num) '.dat'])
                dat=load_rbrs([config.root_dir config.raw_path config.mooring_paths{a} '/rbrs/' fname.name]);
            case 4 % Chipod
                dat.sensor='chipod_no_data_yet'
        end
        save([config.processed_dir config.mooring_paths{a} '/' dat.sensor '_' num2str(serial_num)],'dat');
        else
            % load the matlab files instead
            serial_num=mooring{a}.config.sensornums(b)
                load([config.processed_dir config.mooring_paths{a} '/' config.sensor_types{mooring{a}.config.sensor_type(b)} '_' num2str(serial_num)]);

                % if we're dealing with chipod, but have no data:
                
                if findstr(dat.sensor,'chipod_no_data_yet')
                    dat.time=config.time;
                    dat.T=config.time*NaN;
                end
                 
                % let's filter the data and create a timeseries of T and of
                % var(T)
                delt=mode(diff(dat.time));
                filter_length=config.delt/delt;
                the_filt=hanning(round(filter_length*2));the_filt=the_filt/sum(the_filt);
                % identify periods of no data using diff(dat.time):
                bad_inds=find(diff(dat.time)>1/60/24);
                dat.T(bad_inds)=NaN;
                dat.T(bad_inds+1)=NaN;
                bad_inds2=find(diff(dat.time)<=0);
                if ~isempty(bad_inds2) & findstr(dat.sensor,'sbe56')
                   for abc=1:length(bad_inds2)
                       dat.time=dat.time([ 1:(bad_inds2(abc)-1) (bad_inds2(abc)+1):end]);
                       dat.T=dat.T([ 1:(bad_inds2(abc)-1) (bad_inds2(abc)+1):end]);
                   end
                end
                if mooring{a}.config.sensornums(b)==321
                    dat.T(find(dat.time>datenum('07-Jul-2011')))=NaN;
                end
                
                if mooring{a}.config.sensornums(b)==7818
                    dat.P(dat.P>1416 | dat.P<1410)=1412; % we'll just hardwire it for now.
                    disp('FIX SENSOR 7818 in the future')
                end
                % Byungho said just delete sensor number 39 from the sensor definition
                % file instead....  (at least for the time being...)
                %               if mooring{a}.config.sensornums(b)==39
                %                   dat.T=dat.T*NaN;
                %               end

                % now bin average and compute highpass variance
                T_lowpass=conv2(dat.T,the_filt,'same');
                T_highpass=dat.T-T_lowpass;
                T_var=conv2(T_highpass.^2,the_filt,'same');
                doprint=0
                if doprint
                    colors='rgbm'
                    col=colors(mooring{a}.config.sensor_type(b))
                    plot(dat.time(1:1:end),dat.T(1:1:end),col),hold on
                    plot(dat.time(1:10:end),T_lowpass(1:10:end),'k'),hold on
                    plot(dat.time(1:10:end),100*T_var(1:10:end),'k'),hold on
                end
                % now assign these data to variables in mooring{a}
                
                mooring{a}.T(b,:)=interp1(dat.time,T_lowpass,mooring{a}.time);
                mooring{a}.Tvar(b,:)=interp1(dat.time,T_lowpass,mooring{a}.time);
                if findstr('sbe37',dat.sensor) %dat.config.P % add P if w have it...
                    mooring{a}.P(b,:)=interp1(dat.time,conv2(dat.P,the_filt,'same'),mooring{a}.time);
                    mooring{a}.C(b,:)=interp1(dat.time,conv2(dat.C,the_filt,'same'),mooring{a}.time);
%                    clf,plot(mooring{a}.time,mooring{a}.P), kdatetick
%                    paus
                end
                ylim([0 30])
                kdatetick                
                pause(.01)
        end
      
    end

    xm=mooring{a}
    
    % Now assign pressures to each of the sensors
    base_t=datenum('dec-30-2010');
    
%     if a<5 % this is a T-chain; I think the 
%     sensor_fractions=multcol(ones(n_sensors,len_t),mooring{a}.config.distance_from_top_sbe37/mooring{a}.config.distance_from_top_sbe37(end));
%     delp=mooring{a}.P(end,:)-mooring{a}.P(1,:);
%     mooring{a}.P=addcol(multcol(sensor_fractions,delp),mooring{a}.P(1,:));
%     else % this is the stablemoor
       % hardwire the calculation: 
       % now we'll loop through each of the regions where we have senseors bounded by P 

       
       if a==5 % got some problems here to take care of for the stablemoor first....
         for c=1:2
               microcat_inds=find(mooring{a}.config.sensor_type==2)
               % fix the time period when the lower P sensors went wacko:
               if c==1
                   tinds=find(mooring{a}.time>734716);
                   these_inds=microcat_inds(end-2):microcat_inds(end);
               else
                   tinds=find(mooring{a}.time>734716);
                   these_inds=microcat_inds(end-4):microcat_inds(end-2);
               end
                   
           n_sens=length(these_inds);
           dist_down_wire=mooring{a}.config.distance_from_top_sbe37(these_inds)-mooring{a}.config.distance_from_top_sbe37(these_inds(1));
           sensor_fractions=multcol(ones(n_sens,len_t),dist_down_wire/(mooring{a}.config.distance_from_top_sbe37(these_inds(end))-mooring{a}.config.distance_from_top_sbe37(these_inds(1))));
           delp=mooring{a}.P(these_inds(end),:)-mooring{a}.P(these_inds(1),:);
           mooring{a}.P(these_inds,tinds)=addcol(multcol(sensor_fractions(:,tinds),delp(:,tinds)),mooring{a}.P(these_inds(1),tinds));
           end
           tmp=mooring{a}.T(38,:);tmp(tmp>8)=NaN;mooring{a}.T(38,:)=tmp;
           tmp=mooring{a}.T(29,:);tmp(mooring{a}.time>(base_t+214))=NaN;mooring{a}.T(29,:)=tmp;
           tmp=mooring{a}.T(43,:);tmp(mooring{a}.time>(base_t+216))=NaN;mooring{a}.T(43,:)=tmp;
       end
       
       
       % all of the moorings should be ok with the following (but I haven't
       % verified this yet....
       microcat_inds=find(mooring{a}.config.sensor_type==2)
       
       for c=1:(length(microcat_inds)-1)
           these_inds=microcat_inds(c):microcat_inds(c+1);
           n_sens=length(these_inds);
           dist_down_wire=mooring{a}.config.distance_from_top_sbe37(these_inds)-mooring{a}.config.distance_from_top_sbe37(these_inds(1));
           sensor_fractions=multcol(ones(n_sens,len_t),dist_down_wire/(mooring{a}.config.distance_from_top_sbe37(these_inds(end))-mooring{a}.config.distance_from_top_sbe37(these_inds(1))));
           delp=mooring{a}.P(these_inds(end),:)-mooring{a}.P(these_inds(1),:);
           mooring{a}.P(these_inds,:)=addcol(multcol(sensor_fractions,delp),mooring{a}.P(these_inds(1),:));
       end
       
    
       
       
%    end
        
    
mooring{a}.T_is_interpolated=logical(zeros(size(mooring{a}.T)));
mooring{a}.T_is_interpolated(isnan(mooring{a}.T))=1;
mooring{a}.T=interp_missing_data(mooring{a}.T);

        
    figure(a)
    clf
ax(1)=subplot(211)
pcolor(datenum2yday(mooring{a}.time),mooring{a}.P,mooring{a}.T),shading interp,axis ij
ylim([0 1400])
ylims=ylim;caxis([2 25]);
if a<5;ylims(1)=500;caxis([2 10]);end
ylim(ylims), 
colormap jet
jtext('Temperature [C]')
freezeColors
h=my_colorbar
axes(h)
freezeColors

do_second=0
if do_second
subplot(212)
n_wind=(1/12)/config.delt; % let's do a two hour filter
the_wind=hanning(n_wind)'/sum(hanning(n_wind));
mooring{a}.T_2hour=conv2(mooring{a}.T,the_wind,'same');
mooring{a}.dTdz=abs(gradient(mooring{a}.T_2hour',mooring{a}.config.distance_from_top_sbe37,mooring{a}.time)');
displacements=multcol((mooring{a}.Tvar),1./(nanmean(mooring{a}.dTdz')));
pcolor(mooring{a}.time-base_t,mooring{a}.P,real(log10(displacements))),axis ij,shading flat
%caxis([2.7 3.5])
ylim(ylims)

end

eval(['mooring' num2str(a) '=mooring{a};']);
save(['mooring' num2str(a) ],['mooring' num2str(a)],'config')
end
delete(h)

save all_moorings mooring config

%%
names={'/Users/nash/analysis/iwise-11/data/raw/A1/adcps/SN1426/data_mat/SN1426_iwise11_ALL.mat',...
'/Users/nash/analysis/iwise-11/data/raw/A1/adcps/SN14167/data_mat/SN14167_iwise11_ALL.mat',...
'/Users/nash/analysis/iwise-11/data/raw/A1/adcps/SN8122/data_mat/SN8122_iwise11_ALL.mat'}
clf
filts=[2 12]
for cc=1:2

subplot(2,1,cc)
for aba=1:3
    
    load(names{aba})
    ginds=find(Vel.yday>162 & Vel.yday<168);
    
    pcolor(Vel.yday(ginds),Vel.z(1:end-1),conv2(diff(Vel.u(:,ginds)),ones(1,filts(cc))/1/filts(cc),'same')),shading flat,axis ij, hold on
    axis tight
%     subplot(212)
%     pcolor(Vel.yday,Vel.z,Vel.v),shading flat,axis ij, hold on
%     axis tight

    colormap redblue4
    
end
caxis([-.25 .25])
jtext('zonal velocity u [m/s]; needs cleaning')
    my_colorbar
   
end

%linkaxes(ax)
%set(gcf,'renderer','zbuffer')

%dual_print_pdf('T-example','same')

%% let's throw data on a uniform grid:

new_grid=[0:10:1400]';
gridded.z=new_grid;
gridded.time=mooring{a}.time;
gridded.T=NaN*ones(length(gridded.z),length(gridded.time));

for bb=1:length(gridded.time)
    gridded.T(:,bb)=interp1(mooring{a}.P(:,bb),mooring{a}.T(:,bb),gridded.z);
end


%%

subplot(212)
hold on

[c,h]=contour(datenum2yday(gridded.time),gridded.z,gridded.T,[2 3 4 5 6 8 10 12 16 20 24]);
set(h,'color','k')


dual_print_pdf('T-example3','same')


%%
%plot the T-chain data:



for a=1:4
   load(['mooring' num2str(a)])
   mooring{a}=eval(['mooring' num2str(a)])
    figure(10+a)
    max_p=ceil(max(max(mooring{a}.P))/10)*10;
    
    new_grid=[0:10:1500]+max_p-1500;
    gridded.z=new_grid;
    gridded.time=mooring{a}.time;
    gridded.T=NaN*ones(length(gridded.z),length(gridded.time));

    for bb=1:length(gridded.time)
        if any(~isnan(mooring{a}.P(:,bb))) & mean(diff(mooring{a}.P(:,bb)))~=0
        gridded.T(:,bb)=interp1(mooring{a}.P(:,bb),mooring{a}.T(:,bb),gridded.z);
        end
    end

    
 % do the plotting....
 clf
 pcolor(gridded.time,gridded.z,gridded.T)
 jtext(['Mooring T' num2str(a) ', Temperature [C]  '],.02,1.03)
 caxis([2.25 4.25])
shading interp
axis ij
axis auto
%xlim([datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')])
xlim([datenum('july-18-2011 20:00:00') datenum('july-21-2011 16:00')])
xlim([datenum('july-19-2011 6:00:00') datenum('july-20-2011 8:00')])
kdatetick
set(gcf,'renderer','zbuffer')
my_colorbar
hold on
[c,h]=contour((gridded.time),gridded.z,gridded.T,[2 2.5 3 3.5 4 4.5 5]);
set(h,'color','k')
ylabel('depth [m]')    
ylim([1000 max_p])
print('-r200','-dpng',['pdf/mooring' num2str(a) '_example2.png'])


end



%%
%load and regrid the T-chain data:
warning('off','MATLAB:interp1:NaNinY')


for a=[1:4]
   load(['mooring' num2str(a)])
   mooring{a}=eval(['mooring' num2str(a)])

    eval(['load mooring' num2str(a) '_eps;']);
    eval(['mooring{a}.eps=moor_eps' num2str(a) '.eps;']);
   
    max_p=ceil(max(max(mooring{a}.P))/10)*10;
    
    new_grid=[0:10:1500]+max_p-1500;
    gridded.z=new_grid;
    gridded.time=mooring{a}.time;
    gridded.T=NaN*ones(length(gridded.z),length(gridded.time));
    gridded.eps=NaN*ones(length(gridded.z),length(gridded.time));
    gridded.max_p=max_p;
    
    for bb=1:length(gridded.time)
        if any(~isnan(mooring{a}.P(:,bb))) & mean(diff(mooring{a}.P(:,bb)))~=0
        gridded.T(:,bb)=interp1(mooring{a}.P(:,bb),mooring{a}.T(:,bb),gridded.z);
        gridded.eps(:,bb)=interp1(mooring{a}.P(:,bb),mooring{a}.eps(:,bb),gridded.z,'nearest');
        end
    end

grd{a}=gridded;
    
    
end


save all_gridded grd
%%

%xlim([datenum('july-18-2011 20:00:00') datenum('july-21-2011 16:00')])
%xlim([datenum('july-19-2011 6:00:00') datenum('july-20-2011 8:00')])
% tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];


tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];
tlims=[datenum('june-15-2011 16:00:00') datenum('june-16-2011 20:00')];

tlims=[datenum('june-14-2011 16:00:00') datenum('june-20-2011 20:00')];
tlims=[datenum('june-14-2011 16:00:00') datenum('june-18-2011 20:00')];
% tlims=[datenum('july-6-2011 16:00:00') datenum('july-10-2011 20:00')];
% tlims=[datenum('june-16-2011 02:00:00') datenum('june-16-2011 14:30')];
tlims2=[datenum('june-15-2011 12:00:00') datenum('june-16-2011 14:00')];
%tlims2=[datenum('june-14-2011 16:00:00') datenum('june-18-2011 20:00')];
make_eg_plot=1;
for a=3 %1:4
   figure(10+a)
   gridded=grd{a};
   max_p=gridded.max_p-20
   ginds=find(gridded.time>tlims(1) & gridded.time<tlims(2));
 % do the plotting....
 clf
    hhh(1)=subplot(211)
 pcolor(gridded.time(ginds),gridded.z,gridded.T(:,ginds))
if make_eg_plot
    jtext(['(b) Temperature [C]  '],.02,1.05)
else
    jtext(['Mooring T' num2str(a) ', Temperature [C]  '],.02,1.05)
end
caxis([2.25 4.25])
shading interp
axis ij
axis auto
xlim(tlims2)
kdatetick
%datetick('keeplimits')
set(gca,'linewidth',1)
set(gcf,'renderer','zbuffer')
hc1=my_colorbar([],[],'T [^oC]')
hold on
[c,h]=contour(gridded.time(ginds),gridded.z,gridded.T(:,ginds),[2 2.5 2.7 3 3.5 4 4.5 5]);
set(h,'color','k')
ylabel('depth [m]')    
ylim([1000 max_p])


hhh(2)=subplot(212)

tmp= gridded.eps(:,ginds);
tmp(tmp<1e-10)=NaN;
pcolor(gridded.time(ginds),gridded.z,log10(tmp))
set(gca,'linewidth',1)
if make_eg_plot
    jtext(['(c) log_{10} \epsilon [W/kg]'],.02,1.03)
else
    jtext(['Mooring T' num2str(a) ', log_{10} \epsilon [W/kg]'],.02,1.03)
end
 caxis([-9 -3.5])
shading flat
axis ij
axis auto
xlim(tlims2)
kdatetick
%datetick('keeplimits')
set(gcf,'renderer','zbuffer')
hc2=my_colorbar([],[],'log_{10}\epsilon [W/kg]')
hold on
[c,h]=contour(gridded.time(ginds),gridded.z,gridded.T(:,ginds),[2 2.5 2.7 3 3.5 4 4.5 5]);
set(h,'color','k')
ylabel('depth [m]')    
ylim([1000 max_p])

%print('-r200','-dpng',['pdf/mooring' num2str(a) '_example2.png'])
end
% add BT tide
axes(hhh(2))
resize_gca([0 -.03 -.04])
resize_gca([-.04 -.03],hc2)
axes(hhh(1))
xlabel('')
resize_gca([0 -.10 -.04])
resize_gca([-.04 -.10],hc1)
hpos=get(gca,'position')
axes('position',[hpos(1) hpos(2)+hpos(4)+.06 hpos(3) .08])
%subplot(6,1,1,'align')
gd_inds=find(btvel.time>tlims(1) & btvel.time<tlims(2));
ha=area(btvel.time(gd_inds),btvel.U_tides(gd_inds))
set(gca,'linewidth',1)
set(ha,'linewidth',2,'facecolor',[.7 .7 1])
%hold off, tmp=plot(btvel.time(gd_inds),0*gd_inds,'k'),hold on
%ha=plot(btvel.time(gd_inds),btvel.U_tides(gd_inds))
%set(ha,'linewidth',2,'color',[.7 .7 1])
%%%ha=plot(btvel.time,btvel.U_tides)
%%%set(ha,'linewidth',2)
axis tight
xlim(tlims2)
kdatetick
%datetick('keeplimits')
ylabel('{\itU_{bt}} [m/s] ')
%ylabel('{U_{bt}} [m/s] ')
xlabel('')
jtext(['(a) barotropic velocity [m/s]'],.02,1.19)
%paus
if make_eg_plot
 dual_print_pdf2('TOS_example_figure','same')
 %   print -r200 -dpng eps/TOS_example_figure.png
end

% 
%dual_print_pdf2('xx','same')

%%
% First let's compute time-averages and minimum detectable dissipation rates.
figure(201),clf
cols='rgbm';
for a=1:4
   x=grd{a};
   tmp=x.eps;
   x.eps_mean=nanmean(x.eps,2);
   tmp(tmp<2e-11)=NaN;
   x.eps_mean_estimates=nanmean(tmp,2);
   
   x.eps_noise=nanmin(tmp,[],2);
   
   grd{a}=x;
    
   h(1)=plot(log10(x.eps_noise),x.z,[cols(a) '--'])
   hold on
   h(2)=plot(log10(x.eps_mean),x.z,[cols(a)],'linewidth',1.5)
   h(3)=plot(log10(x.eps_mean_estimates),x.z,[cols(a) '--'],'linewidth',1.5)
   axis ij
   
end
title('red=T1, green=T2, blue=T3, magenta=T4')
legend(h,{'detection limit - min eps','deployment-mean','mean of non-zero overturns'})
xlim([-8 -4])
ylim([600 2300])
xlabel('log_{10} \epsilon [W/kg]')
dual_print_pdf2('thorpe_detection_limits','same')

%%

%load('byungho_Velhpass40h.mat')
clf
x=Velhpass40h;
x.samplerate=1./mean(diff(x.time));
di_filt={'l1.2','h0.8'};
semi_filt={'l2.5','h1.6'};

di_filt={'l1.4','h0.8'};
semi_filt={'l5','h1.4'};

x.U_d1=filter_series(x.Ubt,x.samplerate,semi_filt,4)
x.U_d2=filter_series(x.Ubt,x.samplerate,di_filt,4)
x.U_tides=x.U_d1+x.U_d2;
x.V_d1=filter_series(x.Vbt,x.samplerate,di_filt,4)
x.V_d2=filter_series(x.Vbt,x.samplerate,semi_filt,4)
x.V_tides=x.V_d1+x.V_d2;

x.KE=.5*(x.U_tides.^2+x.V_tides.^2)
x.KE_lp=filter_series(x.KE,x.samplerate,{'l.3'},4)


subplot(211)
plot(x.time,x.Ubt,x.time,x.U_d1,x.time,x.U_d2,x.time,x.U_tides), hold on
kdatetick
subplot(212)
plot(x.time,x.Ubt.^2,x.time,x.U_d1.^2,x.time,x.U_d2.^2,x.time,x.U_tides.^2), hold on
plot(x.time,5*x.KE_lp,'k','linewidth',2)
btvel=rmfield(x,{'Ubc','Vbc'});

save barotropic_vels btvel

%% Now let's look at the temporal cycle of dissipation at all stations.

% Let's compute daily averages, including those from Andy.

new_times=[datenum('14-Jun-2011'):datenum('15-Aug-2011')];
% these correspond roughly to MP times.
del_t=mean(diff(new_times));

% first, the bt_KE is:
load barotropic_vels
KE=interp1(btvel.time,btvel.KE_lp,new_times);

load MP_sn103_overturns_dens_grid
grd{5}=MPeps;
grd{5}.avg.z=MPeps.z';
grd{5}.avg.time=yearday2date(MPeps.ydayavg,2011);
grd{5}.avg.eps=MPeps.epsavg;
tmp=grd{5}.avg.eps;,tmp(isnan(tmp))=1e-11;
grd{5}.avg.inteps=trapz(grd{5}.avg.z,tmp)*1000;

for b=1:4

    x=grd{b};
    
    init=NaN*zeros(length(new_times),length(x.z))';
    x.avg.time=new_times;
    x.avg.z=x.z;
    x.avg.eps=init;
    for a=1:length(new_times);
        inds=find(x.time>new_times(a)-del_t/2 & x.time<new_times(a)+del_t/2);
        x.avg.eps(:,a)=nanmean(x.eps(:,inds),2);
    end
    x.avg.eps(isnan(x.avg.eps))=1e-11;
    x.avg.inteps=trapz(x.avg.z,x.avg.eps)*1000;
    grd{b}=x;
end

%% Plot a scatterplot of integrated dissipation rates. 
stations={'T1','T2','T3','T4','N2'}
figure(104), clf
set(gcf,'position',[335   821   331   290])
cols='rgbmk'
clear h
the_fac=[8e1 8e1 2e2 5e1 5e2]
the_fac=[1 1  1 1 3]
for b=1:5

    subplot(1,1,1)
   h(b)=loglog((2*KE).^(3/2),grd{b}.avg.inteps/the_fac(b),[cols(b) 'o'],'linewidth',1,'markerfacecolor',.75*[1 1 1 ])%cols(b))
   hold on
%xlim(10.^[-2.3 -1])
xlim(10.^[-3.5 -1.5])
ylim([5e-3 10])
grid on
end  
plot(xlim,100*xlim,'k')
set(gca,'linewidth',1)
legend(h,stations,'location','southeast')
xlabel('{\itU_{bt}^3} [m^2/s^2]')
%ylabel('Daily-averaged integrated dissipation \int\epsilon dz [W/m^2]')
ylabel('\int\epsilon dz, 24-h averaged [W/m^2]')
jtext('\epsilon\sim{\itu_{bt}^3}',.82,.62)

dual_print_pdf2('diss_vs_bt3','same')

%% Plot a scatterplot of integrated dissipation rates. 
stations={'T1','T2','T3','T4','N2'}
figure(104), clf
set(gcf,'position',[335   821   338  305])
cols='rgbmk'
clear h
the_fac=[8e1 8e1 2e2 5e1 5e2]
the_fac=[1 1  1 1 3]
leni=5;
tmp=zeros(size(grd{b}.avg.inteps));
for b=1:leni
tmp=tmp+grd{b}.avg.inteps/the_fac(b)/leni
end  

    subplot(1,1,1)
   h(b)=loglog((2*KE).^(3/2),tmp,[cols(b) 'o'],'linewidth',1,'markerfacecolor',.75*[1 1 1 ])%cols(b))
   hold on
%xlim(10.^[-2.3 -1])
xlim(10.^[-3.5 -1.5])
ylim([1e-2 5])
grid on
plot(xlim,120*xlim,'k')
set(gca,'linewidth',1,'fontsize',14)
xlabel('{\itU_{bt}^3} [m^2/s^2]')
ylabel('\int\epsilon dz, 24-h averaged [W/m^2]')
jtext('\epsilon\sim{\itu_{bt}^3}',.82,.76)

dual_print_pdf2('diss_vs_bt4','same')


%% Plot a scatterplot of integrated dissipation rates. 
stations={'T1','T2','T3','T4','N2'}
figure(104), clf
cols='rgbmk'
clear h
the_fac=[8e1 8e1 2e2 5e1 5e2]

for b=1:5
    subplot(3,2,b)
   h(b)=loglog((2*KE).^(3/2),grd{b}.avg.inteps,[cols(b) 'o'],'markerfacecolor',cols(b))
   hold on
xlim(10.^[-3.5 -1.5])
ylim([5e-3 20])
plot(xlim,the_fac(b)*xlim.^(1))
grid on

end   
legend(h,stations,'location','southwest')
xlabel('U_{bt}^3 [m^3/s^3]')
ylabel('Daily-averaged integrated dissipation \int\epsilon dz [W/m^2]')
jtext('\epsilon \sim u_{bt}^3?',.15,.55)

dual_print_pdf2('diss_vs_bt2','same')


%%
%plot a timeseries of daily-averaged dissipations. 
figure(107)
clf
the_inds=[1 2 5 3 4];
for a=1:5
    subplot(6,1,1+a,'align')
    b=the_inds(a);
    pcolor(grd{b}.avg.time,grd{b}.avg.z,log10(grd{b}.avg.eps));
    caxis([-8.5 -5]), axis ij
    yls=ylim;
    xls=xlim;
    ylim(yls(2)+[-1000 0])
shading flat
kdatetick
end
subplot(6,1,1,'align')
plot(btvel.time,btvel.KE_lp,'linewidth',2)
axis tight
xlim(xls)
kdatetick
ylabel('KE [m^2/s^2')

%%
tlims=[datenum('jun-14-2011 0:00:00') datenum('aug-05-2011 0:00')];
%tlims=[datenum('jun-15-2011 0:00:00') datenum('jun-25-2011 0:00')];

figure(109)
 fs=14;
set(gcf,'position',[203   628   709   501])
clf
subplot(5,1,1,'align')
plot(btvel.time,btvel.U_tides,'linewidth',1)
ylim([-.4 .4])
axis tight
xlim(tlims)
kdatetick
ylabel({'{\itU_{BT}} [m/s] '},'fontsize',fs)
xlabel('')
set(gca,'linewidth',1)
jtext('barotropic velocity',.1,0.89)
%resize_gca([0 -.03 0 .03])
subplot(5,1,2,'align')
plot(btvel.time,btvel.KE_lp,'linewidth',2)
axis tight
xlim(tlims)
kdatetick
ylabel({'Barotropic {\itKE}','[m^2/s^2]'},'fontsize',fs)
xlabel('')
ylim([0 .028])
set(gca,'linewidth',1)
jtext('{\itKE_{BT}} [m^2/s^2]',.01,.10)
hplt=subplot(5,1,3,'align')
hold off
hh(1)=plot(grd{1}.avg.time,(conv2(grd{1}.avg.inteps,[1 2 1]/4,'same')),'m--','linewidth',1)
hold on
hh(2)=plot(grd{2}.avg.time,(conv2(grd{2}.avg.inteps,[1 2 1]/4,'same')),'g--','linewidth',1)
hh(3)=plot(grd{3}.avg.time,(conv2(grd{3}.avg.inteps,[1 2 1]/4,'same')),'b--','linewidth',1)
hh(4)=plot(grd{4}.avg.time,(conv2(grd{4}.avg.inteps,[1 2 1]/4,'same')),'r--','linewidth',1)
plot(grd{3}.avg.time,(conv2(.25*(grd{1}.avg.inteps+grd{2}.avg.inteps+grd{3}.avg.inteps+grd{4}.avg.inteps),[1 2 1]/4,'same')),'k','linewidth',2)
resize_gca([0 -.04 0 .04])
hleg=legend(hh,{'T1','T2','T3','T4'},'location','southeast')
resize_gca([.05 -.06],hleg);
set(gca,'linewidth',1)
ylim([.5e-1 3])
set(gca,'yscale','log','ytick',[.1 1])
xlim(tlims)
jtext('\int {\epsilon}dz [w/m^2]',.02,.13)
kdatetick
ylabel({'integrated','\epsilon [W/m^2]'},'fontsize',fs)
xlabel('')
    subplot(3,1,3,'align')
    b=3
    pcolor(grd{b}.avg.time,grd{b}.avg.z,log10(grd{b}.avg.eps));
    caxis([-8.5 -5]), axis ij
    xlim(tlims)
    ylim([900 2080])
shading flat
kdatetick
set(gca,'linewidth',1)
resize_gca([0 0 0 .02])
my_colorbar([],[],'log_{10}\epsilon [W/kg]');%,'fontsize',fs)
set(gca,'layer','top')
ylabel('depth [m]','fontsize',fs)
xlabel('June-August, 2011','fontsize',fs)
jtext('inferred dissipation rate at mooring T3 [W/kg]',.02,1.07)
set(gcf,'renderer','zbuffer')

jet2=jet(64)
jet2(1:8,:)=multcol(jet2(1:8,:),linspace(.25,1,8));
colormap(jet2)
dual_print_pdf2('epsilon_evolution','same')


%
% now compare to Maarten's model 

%x=load('/Users/nash/analysis/iwise-11/maarten/summaries/dissip_WR_2011b.mat')
load('/Users/nash/analysis/iwise-11/maarten/summaries/WR_dissip.mat')

axes(hplt)
inds=[1 2 4 5];
cols='mg br';
the_filt=round((1./mean(diff(pnt(1).tvec)))*24.8/24);
the_wind=hanning(the_filt*1.5)'./sum(hanning(the_filt*1.5));
big_tmp=zeros(size(pnt(1).EPSI))
for b=1:length(inds)
    a=inds(b);
    tmp=real((pnt(a).EPSI))
    tmp=conv2(tmp,the_wind,'same')
    hold on
plot(pnt(a).tvec(the_filt/1:the_filt/2:end),tmp(the_filt/1:the_filt/2:end),[cols(a) '-.'],'linewidth',1)
hold on
big_tmp=big_tmp+tmp./length(inds);

tocomp{b}.time=round(min(pnt(a).tvec)):round(max(pnt(a).tvec));
tocomp{b}.eps_real=interp1(grd{b}.avg.time,(conv2(grd{b}.avg.inteps,[1 2 1]/4,'same')),tocomp{b}.time);
tocomp{b}.eps_model=interp1(pnt(a).tvec,tmp,tocomp{b}.time);
end
plot(pnt(a).tvec(the_filt/1:the_filt/2:end),big_tmp(the_filt/1:the_filt/2:end),['--k'],'linewidth',2)
dual_print_pdf2('epsilon_evolution2','same')

figure(111)
set(gca,'linewidth',1)
clf
cols='mgbr'
for a=1:4
   hhh(a) =loglog(tocomp{a}.eps_real,tocomp{a}.eps_model,[cols(a) '.'])
    hold on
    ginds=find(~isnan(tocomp{a}.eps_real + tocomp{a}.eps_model));
    jtext(['T' num2str(a) ', ' num2str(mean(tocomp{a}.eps_model(ginds))/mean(tocomp{a}.eps_real(ginds)))],.02,.95-.05*a,'color', cols(a))
end
jtext('model bias',.02 ,.95)
plot(xlim,ylim,'k--')
grid on
hleg=legend(hhh,{'T1','T2','T3','T4'},'location','southeast')

xlabel('24-hr average \int \epsilon dz (observations) [W/m^2]')
ylabel('24-hr average \int \epsilon dz (model) [W/m^2]')
dual_print_pdf2('model-obs-comps','same')

%%

% let's also bring in the ADCP at N2.
load('/Users/nash/analysis/iwise-11/m_files/SN8064_iwise11_ALL.mat')

% Now let's do a little phase averaging:
tlims=[datenum('june-14-2011 16:00:00') datenum('july-20-2011 20:00')];
tlims=[datenum('june-11-2011 16:00:00') datenum('aug-4-2011 20:00')];

% let's do a diurnal case and a semidiurnal



for b=1:2

    
    
for cc=1:4
    xx=grd{cc};

    
    %First, let's find the time of maximum eastward flow:
[dd,inds]=findpeaks(btvel.U_tides)
%diurnal_max_inds=inds(dd>0.17 & btvel.time(inds)>xx.time(1))
if b==1
diurnal_max_inds=inds(dd>0.26 & btvel.time(inds)>(xx.time(1)+1));
xtras='';
elseif b==2
diurnal_max_inds=inds(dd<0.26 & dd>.16 & btvel.time(inds)>(xx.time(1)+1))
xtras='_semi';
end
btvel.D1_time=btvel.time(diurnal_max_inds);
btvel.D1_ind=diurnal_max_inds;

    
 fs=14   
clf
set(gcf,'position',[200 400 740 480])
%plot(diff(btvel.time(diurnal_max_inds)))


subplot(5,1,1,'align')
plot(btvel.time,btvel.U_tides,'k','linewidth',1), hold on
plot(btvel.time(diurnal_max_inds),btvel.U_tides(diurnal_max_inds),'r.','markersize',17)
xlim(tlims),kdatetick
ylabel({'barotropic','velocity [m/s]'},'fontsize',fs)
set(gca,'linewidth',1,'layer','top')
% create some phase averages:
len_z=length(xx.z);
len_phase=length(btvel.D1_ind);
len_time=1020;
t_inds=[1:len_time]-len_time/2;
len_time_adcp=len_time*2/5;
adcp_inds=[1:len_time_adcp]-len_time_adcp/2;
tmp1=nan*ones(len_phase,len_z,len_time);
tmp2=nan*ones(len_phase,len_z,len_time);
tmp3=nan*ones(len_phase,length(Vel.z),len_time_adcp);
tmp4=[]
for a=1:len_phase
    the_ind=max(find(xx.time<btvel.D1_time(a)));
    tmp1(a,1:len_z,1:len_time)=xx.T(:,the_ind+t_inds);
    tmp2(a,1:len_z,1:len_time)=xx.eps(:,the_ind+t_inds);
    the_ind2=max(find(Vel.dtnum<btvel.D1_time(a)));
    tmp3(a,1:length(Vel.z),1:len_time_adcp)=Vel.u(:,the_ind2+adcp_inds);
    try
    tmp4(a,1:len_time_adcp)=btvel.U_tides(btvel.D1_ind(a)+adcp_inds);
    catch
    tmp4(a,1:len_time_adcp)=nan;
    end        
end
p_avg.T=squeeze(nanmean(tmp1,1));
p_avg.eps=squeeze(nanmean(tmp2,1));
p_avg.time=t_inds/30;
p_avg.z=xx.z;
p_avg.time_adcp=adcp_inds/30*5/2;
p_avg.z_adcp=Vel.z;
p_avg.U=squeeze(nanmean(tmp3,1));
p_avg.U_A1=squeeze(nanmean(tmp4,1));

max_p=max(xx.z);
xlabel('')
subplot(2,2,3,'align')
 pcolor(p_avg.time,p_avg.z,p_avg.T);
resize_gca([0 0 0 .09])
 caxis([2.25 4.25])
shading interp
axis ij
axis tight
set(gcf,'renderer','zbuffer')
hhh=my_colorbar,resize_gca([-.01 0 -.01],hhh)
hold on
%[c,h]=contour(p_avg.time,p_avg.z,p_avg.T,[2 2.5 3 3.5 4 4.5 5]);
[c,h]=contour(p_avg.time,p_avg.z,p_avg.T,[2 2.5 2.75 3 3.5 4 4.5 5]);
set(h,'color',.05*[ 1 1 1],'linewidth',1)
ylabel('depth [m]','fontsize',fs)    
ylim([900 max_p])
xlabel('Time relative to ebb [h]','fontsize',fs)
jtext('Temperature [C]',.025,1.05,'color','k')
set(gca,'linewidth',1,'layer','top')

subplot(2,2,4,'align')

tmp=p_avg.eps;
z_inds=find(xx.z>900);
tmp=conv2(tmp,hanning(10)'/sum(hanning(10)),'same');
 pcolor(p_avg.time,p_avg.z,log10(tmp));
resize_gca([0 0 0 .09])

caxis([-8 -4.5])
shading flat
axis ij
axis tight
set(gcf,'renderer','zbuffer')
hhh=my_colorbar,resize_gca([-.01 0 -.01],hhh)
hold on
[c,h]=contour(p_avg.time,p_avg.z,p_avg.T,[2 2.5 2.75 3 3.5 4 4.5 5]);
set(h,'color',.95*[ 1 1 1],'linewidth',1)
ylim([900 max_p])
xlabel('Time relative to ebb [h]','fontsize',fs)
jtext('Inferred dissipation rate [W/kg]',.025,1.05,'color','k')
set(gca,'linewidth',1,'layer','top')

subplot(5,2,4)
tmp=conv2(tmp,hanning(20)'/sum(hanning(20)),'same');
p_avg.mean_diss=nanmean(tmp(z_inds,:));
h=area(p_avg.time,log10(p_avg.mean_diss),-9) %,'k','linewidth',2)
set(h,'edgecolor','r','linewidth',2,'facecolor',[1 .75 .75])
ylim([-8 -4.9])
resize_gca([0 -.02])
ylabel({'log_{10}\epsilon','[W/kg]'},'fontsize',fs)
jtext('Inferred dissipation [W/kg]',.025,1.13)
set(gca,'linewidth',1,'layer','top')
subplot(5,2,3)
%h=pcolor(p_avg.time_adcp,p_avg.z_adcp,p_avg.U) %,'k','linewidth',2)
%caxis([-1 1]),shading flat
tmp=conv2(nanmean(p_avg.U(15:20,:)),hanning(10)'/sum(hanning(10)),'same')
h=area(p_avg.time_adcp,tmp,0);
set(h,'edgecolor','b','linewidth',2,'facecolor',[.75 .75 1])
hold on
%h=plot(p_avg.time_adcp,p_avg.U_A1,'k-','linewidth',2);

ylabel({'velocity','[m/s]'},'fontsize',fs)
resize_gca([0 -.02])
jtext('local velocity [m/s]',.025,1.13)
jtext('downslope',.975,.85,'horizontal','right')
jtext('upslope',.975,.15,'horizontal','right')
if b==1
jtext({'local {\itu_{bottom}}'},.2,.75,'color','b')
%jtext({'E ridge {\itU_{bt}}'},.4,.1,'color','k')
ylim([-.5 1])
else 
jtext({'local {\itu_{bottom}}'},.1,.75,'color','b')
%jtext({'East ridge {\itU_{bt}}'},.27,.1,'color','k')
ylim([-.5 1])    
end
set(gca,'linewidth',1,'layer','top')

dual_print_pdf2(['tidal_composite' xtras '_T' num2str(cc)],'same')
pause(.01)
phase_average{b,cc}=p_avg;
end

end


%%
% Now let's generate a movie of this:

% First we need to get the M
% do this later....

the_case=1;

bottom_depth=[1630 1550 1900 2080 2170];
top_depth=900;
big.z=[];
for a=1:length(bottom_depth)
   big.z=[big.z linspace(top_depth,bottom_depth(a),100)']; 
end
the_times=1:length(phase_average{1,1}.time);
big.x=multcol(ones(size(big.z)),1:5);

init=NaN*ones(5,100,length(the_times));
big.T=init;
big.eps=init;
the_tchains=[1 2 4 5];

for b=1:4
    phase_average{the_case,b}.eps2=conv2(phase_average{the_case,b}.eps,hanning(20)'/sum(hanning(20)),'same');
end

for a=the_times
    
    for b=1:length(the_tchains)
        big.T(the_tchains(b),:,a)=interp1(phase_average{the_case,b}.z,phase_average{the_case,b}.T(:,a),big.z(:,the_tchains(b)));
        big.eps(the_tchains(b),:,a)=interp1(phase_average{the_case,b}.z,phase_average{the_case,b}.eps2(:,a),big.z(:,the_tchains(b)));
    end
end
big.T(3,:,:)=0.5*(big.T(4,:,:)+big.T(2,:,:));



%%

tmp=conv2(nanmean(phase_average{the_case,1}.U(15:20,:)),hanning(10)'/sum(hanning(10)),'same')
big.U_local=interp1(phase_average{the_case,1}.time_adcp,tmp,phase_average{the_case,1}.time);
big.time=phase_average{the_case,1}.time;

MakeQTMovie start evolution.mov
n=0
for a=round(135:1.5:885) %;1:1:length(the_times)
clf
set(gcf,'color',[1 1 1])
fs=14
subplot(211)
hold off
  [c,h]=  contourf(big.x(:,the_tchains),big.z(:,the_tchains),squeeze(big.T(the_tchains,:,a))',linspace(2,6,30));
resize_gca([0 -.09 -.05 -.02])
  set(h,'edgecolor','none')
  set(gca,'fontsize',fs,'linewidth',1,'color',.75*[1 1 1])
caxis([2.2 5])
axis ij
ha=arrow([3 2050],[3+2*big.U_local(a)+.1*sign(big.U_local(a)) 2050],'width',2);
my_colorbar([],[],'T [^oC]')
jtext('Temperature ',.025,.13,'fontsize',fs)
if abs(big.time(a))<10
    jtext(['T=' num2str(big.time(a),2) ' h'],.975,1.07,'horizontalal','right')
else
    jtext(['T=' num2str(big.time(a),3) ' h'],.975,1.07,'horizontalal','right')
end
%arrow([3 2000],[3+2*big.U_local(a) 2000]);
%quiver(3,1800,2*big.U_local(a),0,0);
pos=get(gca,'position');
axes('position',[pos(1) pos(2)+pos(4)+.06 pos(3) .10])
ha=area(big.time(50:1020),big.U_local(50:1020),'clipping','off')
set(ha,'edgecolor','k','facecolor',[.75 .75 1],'linewidth',2,'clipping','off')
hold on
plot(big.time(a),big.U_local(a),'bo','linewidth',6,'clipping','off')
set(gca,'xlim',[-15 15],'ylim',[-.3 .9],'ytick',[0 .4 .8],'xticklabel',...
    {'-15','-10','-5','','5','10','15'})
jtext('{\itU_{bottom}} [m/s] ',.025,.63,'fontsize',fs)
ylabel('u [m/s] ')
% text(11,.7,[num2str(big.U_local(a),2) 'm/s'])
xlabel('Time [h]','verticalal','base')
subplot(212)
hold off
[c,h]=   contourf(big.x(:,the_tchains),big.z(:,the_tchains),squeeze(log10(big.eps(the_tchains,:,a)))',[ -11 linspace(-8,-4,30)]);
resize_gca([0 0 -.05 -.02])
set(h,'edgecolor','none')
  set(gca,'fontsize',fs,'linewidth',1,'color',.75*[1 1 1])
 caxis([-8 -5])
axis ij
%hold on
ha=arrow([3 2050],[3+2*big.U_local(a)+.1*sign(big.U_local(a)) 2050],'width',2);
text(2.8,1920,[num2str(big.U_local(a),2) 'm/s'])
my_colorbar([],[],'log_{10} \epsilon [W/kg]')
jtext('Inferred dissipation rate ',.025,0.13,'fontsize',fs)
%hq=quiver(3,1800,2*big.U_local(a),0,0)
%set(hq,'linewidth',2);
xlabel('distance [km]')
ylabel('depth [m]')
set(gcf,'renderer','zbuffer')
pause(.001)
MakeQTMovie addfigure

end
MakeQTMovie finish


%% plot a 4-panel example

%
set(gcf,'position',[ 310   513   607   616])
%xlim([datenum('july-18-2011 20:00:00') datenum('july-21-2011 16:00')])
%xlim([datenum('july-19-2011 6:00:00') datenum('july-20-2011 8:00')])


tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];
tlims=[datenum('june-15-2011 16:00:00') datenum('june-16-2011 20:00')];

tlims=[datenum('june-14-2011 16:00:00') datenum('june-20-2011 20:00')];
tlims=[datenum('june-14-2011 16:00:00') datenum('june-18-2011 20:00')];
tlims=[datenum('july-6-2011 16:00:00') datenum('july-10-2011 20:00')];
tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];
%tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 20:40')];
clf
for a=1:4
   hhh(a)=subplot(5,1,a+1,'align')
   gridded=grd{a};
   max_p=gridded.max_p
   ginds=find(gridded.time>tlims(1) & gridded.time<tlims(2));
   % do the plotting....
   pcolor(gridded.time(ginds),gridded.z,gridded.T(:,ginds))
   caxis([2.25 4.25])
   shading interp
   axis ij
   axis auto
   xlim(tlims)
   kdatetick
   set(gcf,'renderer','zbuffer')
   hold on
   [c,h]=contour(gridded.time(ginds),gridded.z,gridded.T(:,ginds),[2 2.5 3 3.5 4 4.5 5]);
   set(h,'color','k')
   text(tlims(1)+.01, 1000, ['Mooring T' num2str(a) ],'fontweight','bold','color','w')
   ylabel('depth [m]')
   ylim([920 max_p])   
   resize_gca([0 0 -.04])
  xlabel('')
end
% add BT tide
resizes={[.07 -.01],[.11 -.02],[.07 .06],[0.02 .07]}
for a=1:4
resize_gca([0 resizes{a}(1) 0 resizes{a}(2)],hhh(a))
end
hc1=my_colorbar([],[],'T [^oC]')

xlabel('')
hpos=get(hhh(1),'position')
axes('position',[hpos(1) hpos(2)+hpos(4)+.03 hpos(3) .07])
plot(btvel.time,btvel.U_tides,'k','linewidth',2)
hold on
plot(btvel.time,btvel.U_tides*0,'k--','linewidth',1)
axis tight
xlim(tlims)
kdatetick
ylabel('{\itU_{bt}} [m/s]')
xlabel('')
jtext(['barotropic velocity [m/s]'],.02,.84)
dual_print_pdf2('mooring_example_4panel','same')


%% plot a 6-panel example that includes dissipation.

% load barotropic_vels
% load all_gridded

%
set(gcf,'position',[ 310   513   700   640])
%xlim([datenum('july-18-2011 20:00:00') datenum('july-21-2011 16:00')])
%xlim([datenum('july-19-2011 6:00:00') datenum('july-20-2011 8:00')])


tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];
tlims=[datenum('june-15-2011 16:00:00') datenum('june-16-2011 20:00')];

tlims=[datenum('june-14-2011 16:00:00') datenum('june-20-2011 20:00')];
tlims=[datenum('june-14-2011 16:00:00') datenum('june-18-2011 20:00')];
tlims=[datenum('july-6-2011 16:00:00') datenum('july-10-2011 20:00')];
tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 16:00')];
tlims=[datenum('june-15-2011 20:00:00') datenum('june-16-2011 20:40')];
clf
cmap=jet(64);
cmap2=[.9375*[1 1 1]; cmap ];
colormap(cmap2)
cax=[2.25 4.25];
for a=1:4
   hhh(a)=subplot(5,2,2*a+1,'align')
   gridded=grd{a};
   max_p=gridded.max_p
   ginds=find(gridded.time>(tlims(1)-.05) & gridded.time<(tlims(2)+.05));
   % do the plotting....
   tmp=gridded.T(:,ginds);
   tmp(tmp<cax(1))=cax(1);
   pcolor(gridded.time(ginds),gridded.z,tmp)
   caxis([cax(1)-(diff(cax)/64) cax(2)])
   shading interp
   axis ij
   axis auto
   xlim(tlims)
   kdatetick
   set(gcf,'renderer','zbuffer')
   hold on
   [c,h]=contour(gridded.time(ginds),gridded.z,gridded.T(:,ginds),[2 2.5 3 3.5 4 4.5 5]);
   set(h,'color','k')
   text(tlims(1)+.01, 975, ['Mooring T' num2str(a) ],'fontweight','bold','color','w')
if a==2,text(tlims(1)+.00, 760,'(b) Temperature [^oC]','fontweight','bold','color','k'), end
   if a<4,set(gca,'xticklabel',''), else
   ylabel('depth [m]'),end
   ylim([920 max_p])   
   resize_gca([0 0 -.04])
   xlabel('')
   set(gca,'linewidth',1)
   
   hhh2(a)=subplot(5,2,2*a+2,'align')
   % do the plotting of dissipation....
   tmp= gridded.eps(:,ginds);
%   tmp=conv2(tmp,[1 1 1 1 1]/5,'same');
   %   tmp(tmp<1e-10)=1e-10
   tmp(isnan(tmp))=1e-14;
   int_tmp=trapz(gridded.z,tmp);
   int_dissipation(a,1:length(int_tmp))=int_tmp;
   int_time=gridded.time(ginds);
   pcolor(gridded.time(ginds),gridded.z,log10(tmp))
cax2=[-8.5 -4];
   caxis(cax2)
   shading flat
   axis ij
   axis auto
   xlim(tlims)
   if a==2,text(tlims(1)+.00, 820,'(d) log_{10}dissipation [W/kg]','fontweight','bold','color','k'), end
   kdatetick
   set(gcf,'renderer','zbuffer')
   hold on
   [c,h]=contour(gridded.time(ginds),gridded.z,gridded.T(:,ginds),[2 2.5 3 3.5 4 4.5 5]);
   set(h,'color','k')
   %text(tlims(1)+.01, 1000, ['Mooring T' num2str(a) ],'fontweight','bold','color','w')
   %   ylabel('depth [m]')
   ylabel('')
   set(gca,'yticklabel','')
   set(gca,'linewidth',1)
   if a<4,set(gca,'xticklabel',''),end
   ylim([920 max_p])   
   resize_gca([0 0 -.04])
  xlabel('')

end

% add BT tide
resizes={[.07 -.01],[.11 -.02],[.07 .06],[0.02 .07]}
for a=1:4
resize_gca([-.02 resizes{a}(1) 0.05 resizes{a}(2)],hhh(a))
resize_gca([-.04 resizes{a}(1) 0.05 resizes{a}(2)],hhh2(a))
end
%hc1=my_colorbar([],[],'T [^oC]')

load('/Users/nash/analysis/iwise-11/m_files/SN8064_iwise11_ALL.mat')
% this gives you Vel at N2.

xlabel('')
hpos=get(hhh(1),'position')
hhh(5)=axes('position',[hpos(1) hpos(2)+hpos(4)+.02 hpos(3) .07])
vginds2=find(Vel.dtnum>=(tlims(1)) & Vel.dtnum<=(tlims(2)));
tmp=nanmean(Vel.u(:,vginds2));
tmp=conv2(tmp,[1 1]/2,'same');
ha=area(Vel.dtnum(vginds2),tmp,0);
set(ha,'facecolor',[.875 .875 1],'edgecolor',[0 0 .5],'linewidth',1,'clipping','off')
%plot(Vel.dtnum(vginds2),tmp,'k','linewidth',2)
hold on

vginds=find(btvel.time>=(tlims(1)) & btvel.time<=(tlims(2)));
plot(btvel.time(vginds),btvel.U_tides(vginds),'k','linewidth',2)
hold on
plot(btvel.time(vginds),btvel.U_tides(vginds)*0,'k--','linewidth',1)
ylim([-.25 1])
xlim(tlims)
kdatetick
ylabel('{\itU_{bt}} [m/s]')
xlabel('')
set(gca,'xticklabel','','clipping','off')
jtext(['(a) velocity'],.02,1.21)
jtext({'depth-','averaged'},.17,.58,'horizontalalignment','center','fontsize',10)
jtext({'near-bottom'},.7,.83,'color',[0 0 .5],'fontsize',10,'horizontalalignment','center')
   set(gca,'linewidth',1)

hpos=get(hhh2(1),'position')
hhh2(5)=axes('position',[hpos(1) hpos(2)+hpos(4)+.02 hpos(3) .07])
tplt=1020*nanmean(int_dissipation(3:4,:));
tplt=1020*int_dissipation(3,:);
tplt=conv2(tplt,ones(1,5)/5,'same');
tplt(tplt<1e-2)=1e-2;
har=area(int_time,log10(tplt),-2);
set(har,'facecolor',[1 .5 .5])
hold on
plot(int_time,log10(tplt),'color',.5*[1 0 0],'linewidth',1)
xlim(tlims)
ylim([-2 2.5])
kdatetick
ylabel('log_{10}\int\epsilon{\it dz} [W/m^2]    ')
set(gca,'xticklabel','')
   set(gca,'linewidth',1)
xlabel('')
jtext(['(c) depth-integrated turbulent dissipation'],.02,1.21)

dual_print_pdf2('mooring_example_8panel','same')

delete([hhh(1) hhh2(1)])
resize_gca([0 -.09],hhh(5))
resize_gca([0 -.09],hhh2(5))

hpos=get(hhh(2),'position')
hhh1(5)=axes('position',[hpos(1)+hpos(3)-.12 hpos(2)+hpos(4)+.01 .12 .015])
cmp=linspace(cax(1),cax(2),128);
pcolor(cmp,[0 1],[cmp ; cmp]), 
caxis([cax(1)-(diff(cax)/64) cax(2)])
shading flat, set(gca,'ytick',[],'xaxislocation','top','xtick',[2:.5:5])
hold on, plot([1 1]'*[2.5 3 3.5 4],[0 1]','k')
hpos=get(hhh2(2),'position')
hhh2(5)=axes('position',[hpos(1)+hpos(3)-.12 hpos(2)+hpos(4)+.01 .12 .015])
cmp=linspace(cax2(1),cax2(2),128);
pcolor(cmp,[0 1],[cmp ; cmp]), 
shading flat, set(gca,'ytick',[],'xaxislocation','top','xtick',[-10:2:2])

dual_print_pdf2('mooring_example_6panel','same')
print -dpng -r450 eps/mooring_example_6panel_hr.png

%%
tlims=[datenum('june-11-2011 16:00:00') datenum('aug-4-2011 20:00')];

for cc=1:3
tls{1}=[datenum('june-11-2011 16:00:00') datenum('aug-4-2011 20:00')];
tls{2}=[datenum('june-14-2011') datenum('june-20-2011')
    datenum('june-28-2011') datenum('july-5-2011')
    datenum('July-11-2011') datenum('july-18-2011')
    datenum('July-27-2011') datenum('aug-3-2011')]

tls{3}=[datenum('june-20-2011')  datenum('june-28-2011') 
    datenum('july-5-2011')    datenum('July-11-2011') 
    datenum('july-18-2011') datenum('July-27-2011') ]
   figure(4)
   set(gcf,'position',[320 800 400 300])

load MP_sn103_overturns_dens_grid      
load MP_103
clear avgs
% Now make plot that shows the spatial structure of dissipation
avgs.newz=[0:10:2200];
avgs.newz2=big.z;%(:,[1 2 3 4 5]);
avgs.newx2=big.x;%(:,[1 2 3 4 5]);
zinds=1:length(avgs.newz);
to_ind=[1 2 4 5];
for a=[1 2 3 4]
    gridded=grd{a};
    max_p=gridded.max_p;
    lens=size(tls{cc},1);
    ginds=[]
    for b=1:lens
    ginds=[ginds find(gridded.time>tls{cc}(b,1) & gridded.time<tls{cc}(b,2))];
    end
%    ginds=find(gridded.time>tlims(1) & gridded.time<tlims(2));
    avgs.z(to_ind(a),1:151)=gridded.z;
    avgs.epsilon(to_ind(a),1:151)=nanmean(gridded.eps(:,ginds)');
    avgs.eps_newz(to_ind(a),zinds)=interp1(avgs.z(to_ind(a),:),avgs.epsilon(to_ind(a),:),avgs.newz);
    avgs.eps_newz2(to_ind(a),1:100)=interp1(avgs.z(to_ind(a),:),avgs.epsilon(to_ind(a),:),avgs.newz2(:,a));
    avgs.T(to_ind(a),1:151)=nanmean(gridded.T(:,ginds)');
    avgs.T_newz(to_ind(a),zinds)=interp1(avgs.z(to_ind(a),:),avgs.T(to_ind(a),:),avgs.newz);
    avgs.T_newz2(to_ind(a),1:100)=interp1(avgs.z(to_ind(a),:),avgs.T(to_ind(a),:),avgs.newz2(:,a));
end

MPeps.dateavg=yearday2date(MPeps.ydayavg,2011);

lens=size(tls{cc},1);
    ginds=[]
    for b=1:lens
    ginds=[ginds find(MPeps.dateavg>tls{cc}(b,1) & MPeps.dateavg<tls{cc}(b,2))];
    end

avgs.eps_newz2(3,1:100)=interp1(MPeps.z(1:340),nanmean(MPeps.epsavg(1:340,ginds)')/3,avgs.newz2(:,a));

avgs.T_newz2(3,1:100)=NaN; %interp1(MP.z(1:340),nanmean(MP.t(1:340,:)'),avgs.newz2(:,a));    


avgs.eps_toplot=interp_missing_data(avgs.eps_newz2);
avgs.T_toplot=interp_missing_data(avgs.T_newz2);
% this is one way of doing it...
%ginds=[135:885];
%tmp=permute(big.eps(:,:,ginds),[3 1 2]);
%mean_eps=mean(tmp);

clf
[c,h]=contourf(avgs.newx2,avgs.newz2,log10(avgs.eps_toplot'),30),axis ij
set(h,'edgecolor','none')
hold on
ttinds=[1 2 4 5];
[c,h]=contour(avgs.newx2(:,ttinds),avgs.newz2(:,ttinds),log10(avgs.T_toplot(ttinds,:)'),10),axis ij
set(h,'edgecolor','k','linewidth',1)
caxis([-7.5 -5.7])
set(gca,'color',.875*[1 1 1])
text(1,860,'T1','horizontalal','center')
text(2,860,'T2','horizontalal','center')
text(3,860,'N2','horizontalal','center')
text(4,860,'T3','horizontalal','center')
text(5,860,'T4','horizontalal','center')
set(gca,'linewidth',1)
jtext({'time-average','dissipation rate'},.03,.12)
xlabel('cross-ridge distance [km] ')
ylabel('depth [m]')
resize_gca([0 0 -.1])
my_colorbar([],[],'log_{10}\epsilon [W/kg]')
dual_print_pdf2(['time-average-epsilon_' num2str(cc)],'same')
end


