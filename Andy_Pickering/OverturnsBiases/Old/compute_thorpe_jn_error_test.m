addpath /Users/nash/analysis/iwise-11/matlab/ctd/ctd_proc2

% this needs to get run after we have already run load_moorings.
%load all_gridded
for c=[3]

xx=grd{c} 
    
badi=find(isnan(xx.T));
xx.T2=interp_missing_data(extrapolate_data(xx.T,'both'),10);
%xx.T2=extrapolate_data(interp_missing_data(xx.T)','both')';
%xx.T2=(interp_missing_data(xx.T,10));
%paus
xx.S=34.604-.045*(xx.T2-2.5);
xx.eps_from_interped=NaN*xx.S;
h = waitbar(0,'Please wait...');

for a=1:length(xx.time)
    ind=a;
    if mean(xx.T2(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.z',xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
        xx.eps_from_interped(:,ind)=Epsout;
    end
    if mod(a,1000)==0
        waitbar(a/length(xx.time),h)
end
end

%eval(['moor_eps' num2str(c) '.eps=xx.eps;']);
%eval(['save mooring' num2str(c) '_eps moor_eps' num2str(c) ';']);

delete(h)

grd{c}=xx;
end



%%

new_times=[datenum('14-Jun-2011'):datenum('15-Aug-2011')];
del_t=mean(diff(new_times));



for b=3 % 1:4

    x=grd{b};
    
    init=NaN*zeros(length(new_times),length(x.z))';
    x.avg.time=new_times;
    x.avg.z=x.z;
    x.avg.eps2=init;
    for a=1:length(new_times);
        inds=find(x.time>new_times(a)-del_t/2 & x.time<new_times(a)+del_t/2);
        x.avg.eps2(:,a)=nanmean(x.eps_from_interped(:,inds),2);
    end
    x.avg.eps2(isnan(x.avg.eps2))=1e-11;
    x.avg.inteps2=trapz(x.avg.z,x.avg.eps2)*1000;
    grd{b}=x;
end





%%
figure(103)
hl(1)=subplot(211)
imagesc(log10(xx.eps_from_interped))%,shading flat
caxis([-8 -4])
hl(2)=subplot(212)
imagesc(log10(xx.eps))%,shading flat
caxis([-8 -4])
linkaxes(hl)


%%
% 
% for a=1:4
%     
%     eval(['load mooring' num2str(a) ';']);
%     xx=eval(['mooring' num2str(a) ]);
%     eval(['load mooring' num2str(a) '_eps;']);
%     eval(['xx.eps=moor_eps' num2str(a) '.eps;']);
%     moor{a}=xx;
% end
% 
%%
% Now let's test to see if the zig-zagging that the MP does could cause
% some problems with the analysis


addpath /Users/nash/analysis/iwise-11/matlab/ctd/ctd_proc2

n=0
%for tcase=[-1 0 1] 
for tcase=[-1 -.66 -.33 0 .33 .66 1] 
n=n+1;
% this needs to get run after we have already run load_moorings.
%load all_gridded
for c=[3]

xx=grd{c} 
    
x_oset=round(tcase*45*(linspace(-.5,.5,length(xx.z))));
badi=find(isnan(xx.T));
xx.T2=interp_missing_data(extrapolate_data(xx.T,'both'),10);

for a=1:length(xx.z)
    % now let's do some re-assigning:
    maxo=max(abs(x_oset));
    the_inds=(1+maxo):(length(xx.time)-maxo);
    xx.T2(a,the_inds)=xx.T2(a,the_inds+x_oset(a));
end


end
    
%xx.T2=extrapolate_data(interp_missing_data(xx.T)','both')';
%xx.T2=(interp_missing_data(xx.T,10));
%paus
xx.S=34.604-.045*(xx.T2-2.5);
xx.eps_from_interped=NaN*xx.S;
h = waitbar(0,'Please wait...');

for a=10000:27000; %length(xx.time)
    ind=a;
    if mean(xx.T2(:,ind)<6)
        [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns_discrete(xx.z',xx.T2(:,ind),xx.S(:,ind),35.8,0,1,1e-5,0);
        xx.eps_from_interped(:,ind)=Epsout;
    end
    if mod(a,1000)==0
        waitbar(a/length(xx.time),h)
end
end

test_xx{n}=xx;


%eval(['moor_eps' num2str(c) '.eps=xx.eps;']);
%eval(['save mooring' num2str(c) '_eps moor_eps' num2str(c) ';']);

delete(h)

end



save -v7.3 test_xx.mat test_xx



%%
clf

labs={'1500-m zig over 90 min','no zigging or zagging','1500 m zag over 90 min'} 
labs2={'zig','no-zig, no-zag','zag'}
cases=[1 4 7]
for a=1:3
    ginds=[11000:16000];
    hh(a)=subplot(3,1,a,'align')
    imagesc(test_xx{cases(a)}.time(ginds),test_xx{cases(a)}.z,log10((test_xx{cases(a)}.eps_from_interped(:,ginds))))
    resize_gca([0 0 -.35])
    caxis([-8 -4])
    jtext(labs{a})
    kdatetick
    my_colorbar
    hold on
    [c,h]=contour(test_xx{cases(a)}.time(ginds),test_xx{cases(a)}.z,test_xx{cases(a)}.T2(:,ginds),[2 2.5 3 3.5 4 4.5 5]);
    set(h,'color','k')
    
    subplot(3,3,a*3,'align')
    plot(log10(nanmean(test_xx{cases(a)}.eps_from_interped(:,11000:16000)')),-test_xx{cases(a)}.z)
xlim([-8 -5.5])
ylim([-2080 -580])
xlabel('log_{10} \epsilon [W/kg')
   jtext(labs{a})
   jtext(labs2{a},.5, .7,'color','b')

end
    jet2=jet;
    jet2(1,:)=.875*[1 1 1];
    colormap(jet2)
    linkaxes(hh)
    
    
    subplot(3,3,3)
    hold on
    plot(log10(nanmean(test_xx{4}.eps_from_interped(:,11000:16000)')),-test_xx{2}.z,'k')
jtext('no-zig,no-zag',.4, .4)
    subplot(3,3,9)
    hold on
    plot(log10(nanmean(test_xx{4}.eps_from_interped(:,11000:16000)')),-test_xx{2}.z,'k')
jtext('no-zig,no-zag',.4, .4)

   orient landscape
%   print -dpng pdf/epsilon-zig-zagging3.png
   
    

%%
  clf
cols='rgbkbgr';
the_wind=hamming(40*30)'/sum(hamming(40*30));
for a=1:7
   hc(a)= plot(test_xx{a}.time,log10(nanmean(conv2(test_xx{a}.eps_from_interped,the_wind,'same'))),cols(a))
hold on
end
axis tight
kdatetick
legend({'90 min','60 min','30 min','0 min','30 min','60 min','90 min'},'location','southwest')
jtext('depth-averaged \epsilon vs. time to acquire each 1500-m profile')
ylabel('depth-mean \epsilon at T3 [W/kg]')
orient portrait
   print -dpng pdf/epsilon_MP-zigzagging_time-evolution.png


