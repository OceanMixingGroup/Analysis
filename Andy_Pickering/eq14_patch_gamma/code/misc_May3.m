%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% trying to find criteria to exclude convective regions near surface in
% EQ14
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(211);
ezpc(cham.castnumber,cham.P,real(log10(cham.N2)))
caxis([-6 -2])
cb=colorbar
ylim([0 100])
cb.Label.String='log_{10}N^2';cb.FontSize=16;
hline(80,'--')

ax2=subplot(212);
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
caxis([-11 -4])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';cb.FontSize=16;
ylim([0 100])
hline(80,'--')

thresh=0.04
zml=nan*ones(1,length(cham.castnumber));
for i=1:length(cham.castnumber)
%    id=find( cham.SIGMA(:,i)>(cham.SIGMA(1,i)+thresh) );
%if ~isempty(ig)
try
ig = find(~isnan(cham.SIGMA(:,i)));
sigma_surf = cham.SIGMA(ig(1),i);
    id=find( cham.SIGMA(:,i)>(sigma_surf+thresh) );
    if ~isempty(id)
    zml(i)=cham.P(id(1),i);
    end 
end
end

hold on
plot(cham.castnumber,zml,'k')
subplot(211)
hold on
plot(cham.castnumber,zml,'k')


linkaxes([ax1 ax2])

%% save file with depth for each cast to use in thresholding

zml_cnum=cham.castnumber;
save('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat','zml','zml_cnum')


%%