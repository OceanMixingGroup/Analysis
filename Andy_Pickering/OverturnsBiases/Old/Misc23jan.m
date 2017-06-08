
%%


close all

Nm='probability'
%Nm='count'

ht=histogram(Lot_Allcases_True(:),'DisplayStyle','Stair','Normalization',Nm)
hold on
hr=histogram(Lot_Allcases_Resamp(:),'DisplayStyle','Stair','Normalization',Nm)
legend('true','resamp')

%%

ia=find( isnan(Lot_Allcases_True) & ~isnan(Lot_Allcases_Resamp));
ib=find( ~isnan(Lot_Allcases_True) & isnan(Lot_Allcases_Resamp));

%%

close all

h=histogram(Lttot_Allcases_True(:),'DisplayStyle','Stair','Normalization','probability')
hold on
histogram(Lttot_Allcases_Resamp(:),'DisplayStyle','Stair','Normalization','probability')
legend('true','resamp')
xlabel('m')

%%

Eps_Allcases_Resamp(log10(Eps_Allcases_Resamp)<-10)=nan;
Eps_Allcases_True(log10(Eps_Allcases_True)<-10)=nan;

Nm='probability'
%Nm='count'

close all

h=histogram(log10(Eps_Allcases_True(:)),'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(log10(Eps_Allcases_Resamp(:)),'DisplayStyle','Stair','Normalization',Nm)

legend('true','resamp')
xlabel('log_{10} \epsilon','fontsize',16)
grid on
ylabel(Nm,'fontsize',16)
%%

figure(1);clf
histogram(log10(Eps_Allcases_Resamp(:)),'DisplayStyle','Stair','Normalization',Nm)
hold on
histogram(log10(Eps_Allcases_Resamp(ia)),'DisplayStyle','Stair','Normalization',Nm)
%legend('
%%

log10(nanmean(Eps_Allcases_Resamp(:)))
nanmean(log10(Eps_Allcases_Resamp(:)))
log10(nanmean(Eps_Allcases_True(:)))
nanmean(log10(Eps_Allcases_True(:)))
%%