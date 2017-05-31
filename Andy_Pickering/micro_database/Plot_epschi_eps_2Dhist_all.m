%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_epschi_eps_2Dhist_all.m
%
%
%----------------
% 5/25/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all


eps_thresh = 0
eps_floor  = -10

my_data_dir = '/Users/Andy/Google Drive/ChiCalculations/data/'
fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/micro_database/figures'

figure(2);clf
agutwocolumn(1)
wysiwyg

hs=[] ; % collect handles for legend
iax=1
for i = [1:3,5,10]
    
    % close all
    clear tnm K presK Kt hrp hrp96 hrp97
    
    clear mix
    mix = Load_micro_data_mat(i,eps_thresh,eps_floor) ;
    
    
    
    if eps_thresh==1
        %        ib = find(log10(eps)<eps_floor);
        %        eps(ib)=nan;
        %        eps_chi(ib)=nan;
        %        chi(ib)=nan;
        axlims=[eps_floor -6];
    else
        axlims=[-12 -6];
        
    end
    
    
    figure(2)
    subplot(3,2,iax)
    histogram2( real(log10(mix.eps)), real(log10(mix.eps_chi)),'Xbinedges',[-12:0.15:-5],'Ybinedges',[-12:0.15:-5],'DisplayStyle','tile','Normalization','pdf')
    xlim(axlims)
    ylim(axlims)
    xvec = linspace(-12,-5,100);
    hold on
    plot(xvec,xvec,'k--')
    plot(xvec,xvec-1,'r--')
    plot(xvec,xvec+1,'r--')
    title(mix.project)
    xlabel('log_{10}[\epsilon]','fontsize',16)
    ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)
    
    iax = iax + 1;
    
end % i



%~~ add EQ14 data

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')
eps_chi = cham.N2 .* cham.CHI / 2 / 0.2 ./ (cham.DTDZ.^2) ;

if eps_thresh==1
    ib = find(log10(cham.EPSILON)<-8.5);
    %ib = find(cham.P<80);
    cham.EPSILON(ib)=nan;
    eps_chi(ib)=nan;
end


subplot(3,2,6)
histogram2( real(log10(cham.EPSILON)), real(log10(eps_chi)),'Xbinedges',[-12:0.15:-5],'Ybinedges',[-12:0.15:-5],'DisplayStyle','tile','Normalization','pdf')
xvec = linspace(-12,-5,100);
hold on
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
title('EQ14')
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)
xlim(axlims)
ylim(axlims)


% save figure
if eps_thresh==1
    print( fullfile( fig_dir,['epschi_eps_2Dhist_ALL_eps_thresh.png']), '-dpng')
else
    print( fullfile( fig_dir,['epschi_eps_2Dhist_ALL.png']), '-dpng')
end
%%