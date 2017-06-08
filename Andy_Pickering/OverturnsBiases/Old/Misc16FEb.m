%%
%
% Misc16Feb.m
%
% Compare true and resampled quantities
%
%%


clear ; close all
cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

testnum=1
%fname=fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]);
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3_Test' num2str(testnum) '_minOT_50_AllCases'])


d=REsamp.d(:);
dt=REsamp.d_true(:);
figure(1);clf
plot(dt,d,'.')
grid on
gridxy
axis equal
%xlim(500*[-1 1])
%ylim(500*[-1 1])
%%
dd=d-dt;
ig=find(dd~=0);
figure(2);clf
histogram(dd(ig))

%%

L=REsamp.d(:);
Lt=REsamp.d_true(:);
figure(1);clf
plot(Lt,L,'.')
grid on
gridxy
axis equal
%xlim(500*[-1 1])
%ylim(500*[-1 1])
%%
dL=L-Lt;
ig=find(dL~=0);
figure(2);clf
histogram(dL(ig),'DisplayStyle','stair')

nanmean(dL(:))
nanmean(dL(ig))

nanmedian(dL(:))
nanmedian(dL(ig))

%%
N=REsamp.Otnsq_out(:);
Nt=REsamp.Otnsq_out_true(:);
figure(1);clf
plot(Nt,N,'.')
%axis equal
grid on

pe=(N-Nt)./Nt *100;
figure(2);clf
histogram(pe,'DisplayStyle','stair','Normalization','prob')
xlim([-200 800])
grid on
gridxy
%histogram(N-Nt)
%%

n2=REsamp.n2(:);
n2t=REsamp.n2_true(:);
figure(1);clf
%ezpc(REsamp.timeall,REsamp.z,REsamp.n2)
%plot(REsamp.n2_true(:),REsamp.n2(:),'.')
%histogram((n2-n2t)./n2t*100,50)
%histogram(n2-n2t,100)
histogram(n2t,100,'DisplayStyle','stair','Normalization','pdf')
hold on
histogram(n2,100,'DisplayStyle','stair','Normalization','pdf')
%%

clear ; close all
cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases

%fname=fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]);
load '/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3_Test3_minOT_50_AllCases'

%
pm=nan*ones(size(REsamp.eps,1),size(REsamp.eps,3));
pmn2=pm;
for ind=1:size(REsamp.eps,3)
   pm(:,ind)=nanmean(REsamp.eps(:,:,ind),2); 
   pmn2(:,ind)=nanmean(REsamp.n2(:,:,ind),2); 
end

figure(1);clf
semilogx(nanmean(pm,2),REsamp.z)
hold on
semilogx(nanmean(pm(:,1:50),2),REsamp.z)
axis ij

%%

n2=nan*ones(size(REsamp.eps,1),size(REsamp.eps,3));
n2t=n2;
for ind=1:size(REsamp.eps,3)
   n2(:,ind)=nanmean(REsamp.n2(:,:,ind),2); 
   n2t(:,ind)=nanmean(REsamp.n2_true(:,:,ind),2); 
end

figure(1);clf
plot(nanmean(n2t,2),REsamp.z)
hold on
plot(nanmean(n2,2),REsamp.z)
axis ij

%%

figure(1);clf
ezpc(REsamp.timeall(1,:),REsamp.z, log10(nanmean(REsamp.eps,3)))
colorbar
caxis([-9 -3])
%%


clear ; close all
cd /Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases
figure(1);clf

for testnum=[3 1 2 4]
%fname=fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_Case' num2str(whcase) '_minOT_' num2str(minOT)]);
load(['/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases/Data/Tchain3_Test' num2str(testnum) '_minOT_50_AllCases'])

clear N Nt pe
N=REsamp.Otnsq_out(:);
Nt=REsamp.Otnsq_out_true(:);

% plot(Nt,N,'.')
% %axis equal
% grid on

pe=(N-Nt)./Nt *100;
%figure(2);clf
histogram(pe,'DisplayStyle','stair','Normalization','pdf')
hold on
xlim([-200 800])
grid on
%gridxy
%histogram(N-Nt)

end
%%
