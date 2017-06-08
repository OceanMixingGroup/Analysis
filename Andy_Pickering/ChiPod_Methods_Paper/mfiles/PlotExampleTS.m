%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotExampleTS.m
%
% * Makes plot for chipod methods paper *
%
% Modified from: process_ctd_chipod_eq14_AP.m
%
%
%------------------------
% 06/02/16 - A.Pickering - apickering@coas.oregonstate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/ChiPod/EQ14/mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/chipod/
% load paths for data and outputs
Load_chipod_paths_EQ14

% load deployment info
Chipod_Deploy_Info_EQ14

% make list of chipod folders (there is one for each cast)
Fold_list=dir(fullfile(chi_data_path,['*cast*']))
disp(['There are ' num2str(length(Fold_list)) ' folders/casts'])

saveplot=1
%%
% loop through casts
for ifold=2%:length(Fold_list)
    
    close all
    
    %% load CTD data for this cast
    clear castname ctd castnumstr
    castname=Fold_list(ifold).name
    castnumstr=castname(5:6)
    
    % load 24hz data
    load(fullfile(CTD_out_dir_24hz,['CTD' castnumstr '.mat']))
    CTD_24hz=struct();
    CTD_24hz.t1=ctd.temp1(:);
    CTD_24hz.s1=ctd.sal1(:);
    CTD_24hz.p=ctd.pressure(:);
    
    % add time to CTD data
    CTD_24hz.datenum=AddTimetoCTDeq14(ctd);
    CTD_24hz.datenum=CTD_24hz.datenum(:);
    
    % load 1m binned data
    load(fullfile(CTD_out_dir_bin,['CTD' castnumstr '_1m.mat']))
    % dnctd and upctd
    datad_1m=struct();
    datad_1m.t1=dnctd.temp1(:);
    datad_1m.s1=dnctd.sal1(:);
    datad_1m.p=dnctd.pressure(:);
    datad_1m.lat=dnctd.lat;
    datad_1m.lon=dnctd.lon;
    
    datau_1m=struct();
    datau_1m.t1=upctd.temp1(:);
    datau_1m.s1=upctd.sal1(:);
    datau_1m.p=upctd.pressure(:);
    datau_1m.lat=upctd.lat;
    datau_1m.lon=upctd.lon;
    %%
    
    % loop through chipod SNs
    for iSN=1%:2
        
        close all
        
        try
            
            % get name of raw chipod file in this folder
            clear chi_file_list whSN
            whSN=ChiInfo.SNs{iSN}
            this_chi_info=ChiInfo.(whSN);
            chi_file_list=dir(fullfile(chi_data_path,Fold_list(ifold).name,['*cast*.' whSN(3:end)]))
            
            if ~isempty(chi_file_list)
                
                isbig=strcmp('big',ChiInfo.(whSN).InstType)
                
                %~~ specific paths for this chipod
                chi_proc_path_specific=fullfile(chi_proc_path,[whSN]);
                ChkMkDir(chi_proc_path_specific)
                
                chi_fig_path_specific=fullfile(chi_proc_path_specific,'figures')
                ChkMkDir(chi_fig_path_specific)
                %~~
                savedir_cal=fullfile(chi_proc_path_specific,'cal')
                ChkMkDir(savedir_cal)
                
                % do any casts have two chipod files? would need to combine them
                % here
                for ifile=1%:length(chi_file_list)
                    clear fname
                    fname=fullfile(chi_data_path,Fold_list(ifold).name,chi_file_list(ifile).name)
                    
                    % load chipod data (both here are 'big' chipods
                    [data head]=raw_load_chipod(fname);
                    chidat.datenum=data.datenum;
                    len=length(data.datenum);
                    if mod(len,2)
                        len=len-1; % for some reason datenum is odd!
                    end
                    chidat.T1=makelen(data.T1(1:(len/2)),len);
                    chidat.T1P=data.T1P;
                    chidat.T2=makelen(data.T2(1:(len/2)),len);
                    chidat.T2P=data.T2P;
                    chidat.AX=makelen(data.AX(1:(len/2)),len);
                    chidat.AY=makelen(data.AY(1:(len/2)),len);
                    chidat.AZ=makelen(data.AZ(1:(len/2)),len);
                    
                    % plot the chipod data
                    %                     figure;clf
                    %                     agutwocolumn(1)
                    %                     wysiwyg
                    %                     ax = MySubplot(0.1, 0.03, 0.02, 0.06, 0.1, 0.02, 1,4);
                    %
                    %                     axes(ax(1))
                    %                     plot(chidat.datenum,chidat.T1);
                    %                     hold on
                    %                     plot(chidat.datenum,chidat.T2);
                    %                     ylabel('T [V]')
                    %                     datetick('x')
                    %                     grid on
                    %                     legend('T1','T2')
                    %                     title('Raw chipod data (uncalibrated)')
                    %
                    %                     axes(ax(2))
                    %                     plot(chidat.datenum,chidat.T1P);
                    %                     axis tight
                    %                     ylabel('T1P [V]')
                    %                     grid on
                    %                     datetick('x')
                    %
                    %                     axes(ax(3))
                    %                     plot(chidat.datenum,chidat.T2P);
                    %                     axis tight
                    %                     ylabel('T2P [V]')
                    %                     grid on
                    %                     datetick('x')
                    %
                    %                     axes(ax(4))
                    %                     plot(chidat.datenum,chidat.AX,chidat.datenum,chidat.AY,chidat.datenum,chidat.AZ);
                    %                     axis tight
                    %                     ylabel('A [V]')
                    %                     legend('Ax','AY','Az','location','best')
                    %                     datetick('x')
                    %                     grid on
                    %                     xlabel(['Time on ' datestr(floor(nanmin(chidat.datenum))) ])
                    %
                    %                     linkaxes(ax,'x')
                    
                    chidat.castname=castname;
                    
                    % carry over chipod info
                    chidat.Info=this_chi_info;
                    chidat.cal=this_chi_info.cal;
                    az_correction=this_chi_info.az_correction;
                    
                    if strcmp(whSN,'SN1001') && ifold==1
                        % in SN1001, AX is actually AZ for 1st cast onnly?
                        chidat=rmfield(chidat,'AZ');
                        chidat.AZ=chidat.AX;
                    elseif strcmp(whSN,'SN1002')
                        % in SN1002, AY is actually AZ
                        chidat=rmfield(chidat,'AZ');
                        chidat.AZ=chidat.AY;
                    end
                    
                    % align chipod with CTD
                    [CTD_24hz chidat]=AlignChipodCTD(CTD_24hz,chidat,az_correction,1);
                    %                    print('-dpng',fullfile(chi_fig_path_specific,[whSN '_' castname '_Fig2_w_TimeOffset']))
                    
                    % zoom in and plot again
                    xlim([nanmin(chidat.datenum)+1500/86400 nanmin(chidat.datenum)+1600/86400])
                    %                    print('-dpng',fullfile(chi_fig_path_specific,[whSN '_' castname '_Fig3_w_TimeOffset_Zoom']))
                    
                    % Calibrate T and dT/dt
                    [CTD_24hz chidat]=CalibrateChipodCTD(CTD_24hz,chidat,az_correction,1);
                    %                   print('-dpng',fullfile(chi_fig_path_specific,[whSN '_' castname '_Fig4_dTdtSpectraCheck']))
                    
                    % check if T1 calibration is ok
                    clear out2 err pvar
                    out2=interp1(chidat.datenum,chidat.cal.T1,CTD_24hz.datenum);
                    err=out2-CTD_24hz.t1;
                    pvar=100* (1-(nanvar(err)/nanvar(CTD_24hz.t1)) );
                    if pvar<50
                        disp('Warning T calibration not good')
                        fprintf(fileID,' *T1 calibration not good* ');
                    end
                    
                    if this_chi_info.isbig==1
                        % check if T2 calibration is ok
                        clear out2 err pvar
                        out2=interp1(chidat.datenum,chidat.cal.T2,CTD_24hz.datenum);
                        err=out2-CTD_24hz.t1;
                        pvar=100* (1-(nanvar(err)/nanvar(CTD_24hz.t1)) );
                        if pvar<50
                            disp('Warning T2 calibration not good')
                            fprintf(fileID,' *T2 calibration not good* ');
                        end
                    end
                    
                    %~~~~
                    do_timeseries_plot=1;
                    if do_timeseries_plot
                        %%
                        h=ChiPodTimeseriesPlot(CTD_24hz,chidat);
                        axes(h(1))
                        title(['EQ14 - ' castname ' - \chi pod ' whSN  ])%,'interpreter','none')
                        axes(h(end))
                        xlabel(['Time on ' datestr(chidat.datenum(1),'dd-mmm-yyyy')])
                        
                        axes(h(3))
                        ylim([1.5 2.1])
                        ylabel('Acc. [V]')
                        
                        axes(h(5))
                        SubplotLetterMW('Up-looking Sensor',0.3)
                        
                        axes(h(6))
                        SubplotLetterMW('Down-looking Sensor',0.3)
                        
                        linkaxes(h,'x')
                        %%
                        if saveplot==1
                            SetPaperFigPath
                            print('-dpng','-r300',fullfile(figdir,[whSN '_' castname '_Fig5_T_P_dTdz_fspd.png']));
                        end
                        
                    end
                    
                end
            end
        end
    end
end
%%