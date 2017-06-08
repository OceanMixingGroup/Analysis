%%
%
% Misc5Mar.m
%
% Plot histograms for a specific depth range and short time period. See
% what is causes averages to differ.
%
%%

z_range=[1400 1450]

iz=isin(REsamp.z,z_range);

figure(1);clf

ecall=[]
for whc=1:Nshift
    clear idt ec tc zc whdir
    idt=isin(REsamp.timeall(whc,:),xl);

    ec=squeeze(REsamp.eps(:,idt,whc));
    tc=squeeze(REsamp.tsamp(:,idt,whc));
    zc=squeeze(REsamp.zsamp(:,idt,whc));

    ecall=[ecall ec(iz,:)];
end

Nm='pdf'
Nm='count'

histogram(log10(ecall(:)),40,'Normalization',Nm)
hold on
izt=isin(xx2.z,z_range);
histogram(log10(xx2.eps(izt,idtrue)),40,'Normalization',Nm)
legend('sample','true')

etrue=xx2.eps(izt,idtrue);
freqline(log10(nanmean(ecall(:))),'k')
freqline(log10(nanmean(etrue(:))),'r')
xlabel('log_{10}\epsilon','fontsize',16)

%%

Nm='pdf'
%Nm='count'

figure(2);clf
idgs=find(log10(ecall)>-10);
idgt=find(log10(etrue)>-10);
ht=histogram(log10(etrue(idgt)),20,'Normalization',Nm)
hold on
hs=histogram(log10(ecall(idgs)),20,'Normalization',Nm)

legend('true','sampled','location','best')

%freqline(log10(nanmean(ecall(idgs))),'k')
%freqline(log10(nanmean(etrue(idgt))),'r')

freqline(log10(nanmean(ecall(:))),'k')
freqline(log10(nanmean(etrue(:))),'r')

hold on
lims=ylim;
plot(log10(nanmean(ecall(:))),lims(2),'kd','markersize',10)
plot(log10(nanmean(etrue(:))),lims(2),'rs','markersize',10)

xlabel('log_{10}\epsilon','fontsize',16)

%%
