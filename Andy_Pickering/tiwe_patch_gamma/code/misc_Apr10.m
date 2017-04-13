%
%
% Trying to figure out why chi only goes to ~150m for profiles in first
% part of record?
%
% many T profiles seem to be totally constant below ~160m (same depth where
% TP and chi are missing below)?
%
%%


for cnum=20%1:50:4000
    
    try
        
        load(['/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/cham_proc/cal/tw91_' sprintf('%04d',cnum) '_cal.mat'])
        
        figure(1);clf
        set(gcf,'Name',['profile' num2str(cnum)])
        
        ax1=subplot(121)
        plot(cal2.T1,cal2.P,'o-')
        axis ij
        grid on
        
        ax2=subplot(122)
        plot(cal2.TP,cal2.P)
        axis ij
        grid on
        xlim([-10 10])
        
        linkaxes([ax1 ax2],'y')
        pause(1)

    catch
    end % try

end % cnum

%%
