%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CombineResampTestCases_LES.m
%
% Combine all T-chain resampling test cases (entire ensemble) into one structure to
% make it easier to load and analyze.
%
% Modified from CombineResampTestCases.m
%
% Resampled data is made in Resample_LES_OT.m
%
% 20 Feb. 2015 - A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES

savedata=1

dztot=1835;

minOT=50

for  whtest=1:5
    
    clear testnum w_samp
    if whtest==1
        w_samp=0.25; testnum=1 % m/s
    elseif whtest==2
        w_samp=0.5; testnum=2 % m/s
    elseif whtest==3
        w_samp=0.15; testnum=3 % m/s
    elseif whtest==4
        w_samp=0.75; testnum=4 % m/s
    elseif whtest==5
        w_samp=1.0; testnum=5 % m/s
    else
    end
    
    %
    tshift=(2/60/24) % shift each case by a few minutes to create ensemble of all phases
    Tprof=dztot/w_samp/86400; % time for a profile (1 up or down) in days
    Nshift=round(2*Tprof/tshift)% number of shifts to go through all phases
    
    % Load true T-chain data
    clear xx2 REsamp
    %    load(fullfile('Data',['Tchain' num2str(whmoor) '_RecomputedEps_minOT_' num2str(minOT)]))
    
    DatDir=['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/LES/Data/Test' num2str(testnum)]
    
    %~~ first load 1 resampled dataset to get array sizes etc.
    whcase=1
    clear fname x_resamp
    fname=['LES_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
    load(fullfile(DatDir,fname))
    
    % make empty arrays for combined data
    % dimensions are depth X time X whcase
    EmptyMat=nan*ones( length(x_resamp.z),length(x_resamp.time),Nshift);
    REsamp=struct();
    REsamp.eps=EmptyMat;
    REsamp.t=EmptyMat;
    REsamp.Lot=EmptyMat;
    REsamp.Lttot=EmptyMat;
    REsamp.d=EmptyMat;
    REsamp.Otnsq_out=EmptyMat;
    REsamp.n2=EmptyMat;
    
    %    REsamp.start=EmptyMat;
    REsamp.Lt_each=EmptyMat;
    REsamp.Lot_each=EmptyMat;
    REsamp.eps_each=EmptyMat;
    REsamp.Otnsq_each=EmptyMat;
    
    REsamp.zall=nan * ones( length(x_resamp.z) , Nshift ) ;
    REsamp.timeall=nan * ones( Nshift, length(x_resamp.time) ) ;
    
    REsamp.tsamp=EmptyMat;
    REsamp.zsamp=EmptyMat;
    
    % arrays for 'true' data corresponding to resampled profiles
    REsamp.eps_true=EmptyMat;
    REsamp.t_true=EmptyMat;
    REsamp.Lot_true=EmptyMat;
    REsamp.Lttot_true=EmptyMat;
    REsamp.d_true=EmptyMat;
    REsamp.Otnsq_out_true=EmptyMat;
    REsamp.n2_true=EmptyMat;
    REsamp.timeall_true=nan * ones( Nshift, length(x_resamp.time) ) ;
    
    REsamp.Lt_each_true=EmptyMat;
    REsamp.Lot_each_true=EmptyMat;
    REsamp.eps_each_true=EmptyMat;
    REsamp.Otnsq_each_true=EmptyMat;
    
    % get data from each case (different start times/phases)
    hb=waitbar(0,'workin away')
    
    for whcase=1:Nshift
        
        waitbar(whcase/Nshift,hb)
        clear fname x_resamp
        fname=['LES_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]   ;
        load(fullfile(DatDir,fname))
        
        %        REsamp.t        (:,:,whcase)=x_resamp.T;
        
        REsamp.Lot      (:,:,whcase)=x_resamp.Lot;
        REsamp.Lttot    (:,:,whcase)=x_resamp.Lttot;
        REsamp.eps      (:,:,whcase)=x_resamp.eps;
        REsamp.d        (:,:,whcase)=x_resamp.d;
        REsamp.Otnsq_out(:,:,whcase)=x_resamp.Otnsq_out;
        REsamp.n2       (:,:,whcase)=x_resamp.n2;
        REsamp.zall     (:,whcase)  =x_resamp.z;
        REsamp.timeall  (whcase,:)  =x_resamp.time;
        
        %        REsamp.start(:,:,whcase)=x_resamp.start;
        REsamp.Lt_each(:,:,whcase)=x_resamp.Lt_each;
        REsamp.Lot_each(:,:,whcase)=x_resamp.Lot_each;
        REsamp.eps_each(:,:,whcase)=x_resamp.eps_each;
        REsamp.Otnsq_each(:,:,whcase)=x_resamp.Otnsq_each;
        
        REsamp.tsamp(:,:,whcase)=x_resamp.tsampall;
        REsamp.zsamp(:,:,whcase)=x_resamp.zall;
        
        %%~~~~~~~
        %
        %         %~~ Get true profiles corresponding to resampled times (midpoints of
        %         %resampled profiles). This section of code taken from
        %         % PlotResampVsTrue_direct.m.
        %         % find indices of Tchain data corresponding to time of each profile
        %         clear It
        %         It=nan*ones(size(x_resamp.time));
        %         for wht=1:length(x_resamp.time)
        %             clear val
        %             [val,It(wht)]=nanmin(abs(xx2.yday-x_resamp.time(wht))) ;
        %             if val>1/24
        %                 It(wht)=nan;
        %             end
        %         end
        %
        %         REsamp.t_true(:,:,whcase)=xx2.T(:,It);
        %         REsamp.Lot_true(:,:,whcase)=xx2.Lot(:,It);
        %         REsamp.Lttot_true(:,:,whcase)=xx2.Lttot(:,It);
        %         REsamp.eps_true(:,:,whcase)=xx2.eps(:,It);
        %         REsamp.d_true(:,:,whcase)=xx2.d(:,It);
        %         REsamp.Otnsq_out_true(:,:,whcase)=xx2.Otnsq_out(:,It);
        %         REsamp.n2_true(:,:,whcase)=xx2.n2(:,It);
        %         REsamp.zall_true(:,whcase)=xx2.z;
        %         REsamp.timeall_true(whcase,:)=xx2.yday(It);
        %
        %         REsamp.Lt_each_true(:,:,whcase)=xx2.Lt_each(:,It);
        %         REsamp.Lot_each_true(:,:,whcase)=xx2.Lot_each(:,It);
        %         REsamp.eps_each_true(:,:,whcase)=xx2.eps_each(:,It);
        %         REsamp.Otnsq_each_true(:,:,whcase)=xx2.Otnsq_each(:,It);
        %         %
        %
        %         %%~~~~~~~
        
    end % whcase
    
    delete(hb)
    
    % these fields are the same for all cases (start times)
    REsamp.w_samp=x_resamp.w_samp;
    REsamp.z=x_resamp.z;
    REsamp.whmoor='LES';
    REsamp.testnum=testnum;
    REsamp.minOT=minOT;
    REsamp.MakeInfo=['Made ' datestr(now) ' w/ CombineResampTestCases_LES.m in ' version]
    
    if savedata==1
        fname=fullfile(DatDir,['LES_Test' num2str(testnum) '_minOT_' num2str(minOT) '_AllCases'])
        save(fname,'REsamp')
    end
    
    
end % whtest

%%