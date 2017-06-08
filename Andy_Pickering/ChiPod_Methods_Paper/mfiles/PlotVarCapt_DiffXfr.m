%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotVarCapt_DiffXfr.m
%
% * Makes plot for chipod methods paper *
%
% Make a plot showing the % of variance of spectra captured when not using
% a response correction, based on the (historically) measured response functions.
%
% I will do this for 2 different freq ranges: <=7Hz, and <=15Hz. Hopefully
% will show that using only up to 7Hz captures most of variance even with
% unknown xfr functions.
%
%--------------
% 05/27/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Dropbox/transfer_functions/transfer_functions_yq08a.mat')
Nnames=83;
str='yq08a'

names=fields(transfer)
names=names(1:Nnames);

%%

v1=[]
v2=[]
for isens=1:length(names)
    %    isens=10
    whsens=names{isens};
    
    fc=transfer.filter_freq.(whsens);%10;
    n_order=transfer.filter_ord.(whsens);%2;
    freq=transfer.f;%f_obs(1,:);
    spec=ones(size(freq));
    out=1./(1+(freq./fc).^(2*n_order));
    
    figure(1);clf
    hcor=loglog(freq,out,'b-');
    hold on
    loglog(transfer.f,transfer.(whsens),'color',0.7*[1 1 1],'linewidth',1)
    loglog(transfer.smoothed_spec.f,transfer.smoothed_spec.(whsens),'k')
    grid on
    freqline(7)
    freqline(15)
    
    
    id1=isin(freq,[0 7]);
    id2=isin(freq,[0 15]);
    
    %
    df=nanmean(diff(freq));
    var1=sum(ones(size(id1))) / sum(1./out(id1))  * 100;
    var2=sum(ones(size(id2))) / sum(1./out(id2))  * 100;
    
    v1=[v1 var1];
    v2=[v2 var2];
    
    %pause
    
end
%%

figure(1);clf
subplot(121)
h1=histogram(v1)
hold on
subplot(122)
histogram(v2,h1.BinEdges)
%%
figure(1);clf

subplot(121)
boxplot(v1)
ylim([60 100])
title('0<f<10')
ylabel('% variance captured')

subplot(122)
boxplot(v2)
ylim([60 100])
title('0<f<15')

%%
% figdir='/Users/Andy/Cruises_Research/ChiPod/ChiPod_Methods_Paper'
% print(fullfile(figdir,['EQ14_boxplot_PctVar_xfr']),'-dpng')

%% Compute percent of times the variance captured is >90%

id1=find(v1>95);
p1=length(id1)/length(v1)*100

id2=find(v2>95);
p2=length(id2)/length(v2)*100

%%