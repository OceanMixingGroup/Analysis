%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotOldTransFunc.m.
%
%
%
%----------------
% 06/24/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all
saveplot=1
SetParamsPaperPlots
for iyr=3%1:3
    
    switch iyr
        case 1
            load('/Users/Andy/Dropbox/transfer_functions/transfer_functions_yq06b.mat')
            Nnames=20;
            str='yq06b'
        case 2
            load('/Users/Andy/Dropbox/transfer_functions/transfer_functions_yq07a.mat')
            Nnames=35;
            str='yq07a'
        case 3
            load('/Users/Andy/Dropbox/transfer_functions/transfer_functions_yq08a.mat')
            Nnames=83;
            str='yq08a'
    end
    
    names=fields(transfer)
    names=names(1:Nnames);
    
    % Plot all transfer functions, and the one I found works best for EQ14
    % chameleon data.
    figure(1);clf
    agutwocolumn(0.6)
    wysiwyg
    
    for isens=1:length(names)
        
        whsens=names{isens}
        
        fc=transfer.filter_freq.(whsens);%10;
        n_order=transfer.filter_ord.(whsens);%2;
        freq=transfer.f;%f_obs(1,:);
        spec=ones(size(freq));
        out=1./(1+(freq./fc).^(2*n_order));
        hcor=loglog(freq,out,'-','color',0.7*[1 1 1],'linewidth',1);
        hold on
        grid on
        
    end
    
    %     fc=10%transfer.filter_freq.(whsens);%10;
    %     n_order=2%transfer.filter_ord.(whsens);%2;
    %     freq=transfer.f;%f_obs(1,:);
    %     spec=ones(size(freq));
    %     out=1./(1+(freq./fc).^(2*n_order));
    %     hcor=loglog(freq,out,'r--','linewidth',2)
    
    fc=32%transfer.filter_freq.(whsens);%10;
    n_order=2%transfer.filter_ord.(whsens);%2;
    freq=transfer.f;%f_obs(1,:);
    spec=ones(size(freq));
    out=1./(1+(freq./fc).^(2*n_order));
    hdef=loglog(freq,out,'b--','linewidth',1)
    
    
    xlim([1 200])
    ylim([ 10^(-1.5) 10^(0.5)])
    freqline(2);
    freqline(7);
    xlabel('frequency [Hz]','fontsize',16)
    text(2+0.25,10^(0.25),['2Hz'],'fontsize',15)
    text(Params.fmax+0.05,10^(0.25),[num2str(Params.fmax) 'Hz'],'fontsize',15)
    title(['Measured FP07 Transfer Functions '])
    
    if saveplot==1
        SetPaperFigPath
        figname=[str '_TransFunc']
        print(fullfile(figdir,figname) , '-dpng' )
    end
    
end
%%