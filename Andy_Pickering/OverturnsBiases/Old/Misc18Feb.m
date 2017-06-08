%%
%
%
%
%
%%  Make example of how histograms for Lt are skewed when using original data that has same value repeated for each depth in overturn

clear ; close all

testnum=3

load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3_Test' num2str(testnum) '_minOT_50_AllCases'])

figure(1);clf
histogram(REsamp.Lttot_true(:),'Normalization','pdf')
hold on
histogram(REsamp.Lt_each_true(:),'Normalization','pdf')
xlabel('L_T')
title('PDF of Thorpe scales')
%%
Ds='stair'
%Ds='bar'
figure(1);clf
histogram(REsamp.Lt_each_true(:),'Normalization','pdf','DisplayStyle',Ds)
hold on
histogram(REsamp.Lt_each(:),'Normalization','pdf','DisplayStyle',Ds)
legend('true','sampled')
xlabel('L_T')

%%
Ds='stair'
%Ds='bar'
figure(1);clf
histogram(REsamp.Otnsq_each_true(:),'Normalization','pdf','DisplayStyle',Ds)
hold on
histogram(REsamp.Otnsq_each(:),'Normalization','pdf','DisplayStyle',Ds)
legend('true','sampled')
xlim([0 1e-5])
xlabel('N^2')
%%
Ds='stair'
%Ds='bar'
figure(1);clf
histogram(log10(REsamp.eps_each_true(:)),'Normalization','pdf','DisplayStyle',Ds)
hold on
histogram(log10(REsamp.eps_each(:)),'Normalization','pdf','DisplayStyle',Ds)
legend('true','sampled')

%%

ig=find(log10(REsamp.eps_true(:))>-10);
ig2=find(log10(REsamp.eps(:))>-10);

Ds='stair'
%Ds='bar'
figure(1);clf
%histogram(log10(REsamp.eps_true(:)),'Normalization','pdf','DisplayStyle',Ds)
histogram(log10(REsamp.eps_true(ig)),'Normalization','pdf','DisplayStyle',Ds)
hold on
%histogram(log10(REsamp.eps(:)),'Normalization','pdf','DisplayStyle',Ds)
histogram(log10(REsamp.eps(ig2)),'Normalization','pdf','DisplayStyle',Ds)
legend('true','sampled')

%% Try looking at just one depth

Ds='stair'
%Ds='bar'

%[val,iz]=nanmin(abs(REsamp.z-1000))
iz=isin(REsamp.z,[1400 1450])
%iz=isin(REsamp.z,[1000 1200])
REsamp.z(iz) 

e1=REsamp.eps_true(iz,:,:);
e2=REsamp.eps(iz,:,:);

Lt1=REsamp.Lttot_true(iz,:,:);
Lt2=REsamp.Lttot(iz,:,:);

N1=REsamp.Otnsq_out_true(iz,:,:);
N2=REsamp.Otnsq_out(iz,:,:);

ig1=find(log10(e1)>-10);
ig2=find(log10(e2)>-10);
%%
log10(nanmean(e1(:)))
log10(nanmean(e2(:)))

% log10(nanmean(e1(ig1)))
% log10(nanmean(e2(ig2)))
%%
figure(1);clf
agutwocolumn(1)
wysiwyg
ax = MySubplot(0.1, 0.03, 0.075, 0.06, 0.1, 0.05, 2,2);

axes(ax(1))
histogram(log10(e1(ig1)),'Normalization','pdf','DisplayStyle',Ds)
hold on
%Ds='stair'
histogram(log10(e2(ig2)),'Normalization','pdf','DisplayStyle',Ds)
xlabel('log_{10}\epsilon')
legend('true','sampled')

axes(ax(2))
%histogram(Lt1(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(Lt1(ig1),'Normalization','pdf','DisplayStyle',Ds)
hold on
%Ds='stair'
%histogram(Lt2(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(Lt2(ig2),'Normalization','pdf','DisplayStyle',Ds)
xlabel('L_T')
legend('true','sampled')


axes(ax(3))
%histogram(Lt1(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(N1(ig1),'Normalization','pdf','DisplayStyle',Ds)
hold on
%Ds='stair'
%histogram(Lt2(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(N2(ig2),'Normalization','pdf','DisplayStyle',Ds)
xlabel('N^2')
legend('true','sampled')

d1=REsamp.d_true(iz,:,:);
d2=REsamp.d(iz,:,:);



axes(ax(4))
%histogram(Lt1(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(d1(ig1),'Normalization','pdf','DisplayStyle',Ds)
hold on
%Ds='stair'
%histogram(Lt2(:),'Normalization','pdf','DisplayStyle',Ds)
histogram(d2(ig2),'Normalization','pdf','DisplayStyle',Ds)
xlabel('d')
legend('true','sampled')
grid on

%%
figure(1);clf
%plot(Lt1(:),Lt2(:),'.')
loglog(N1(:),N2(:),'.')
%loglog(e1(:),e2(:),'.')
%axis equal
grid on
%%
