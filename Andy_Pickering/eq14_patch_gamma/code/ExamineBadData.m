%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ExamineBadData.m
%
%
%
%-----------------
% 4/26/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7  ;
Params.z_smooth=10;
dz = 2 ;
cnums_to_get = get_cham_cnums_eq14;

%bad_prof=[2282 2283 2391 2762 2953] % bad profiles where temp. was wacky
%cnums_to_get = setdiff(cnums_to_get,bad_prof);

eq14_patches_paths
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200);    

%% Plot N2 vs dTdz, with dots colored by epsilon?

id = find(chipod.P>80);

X = chipod.N2(id,:); X=X(:);
Y = chipod.Tz(id,:); Y=Y(:);
C = chipod.eps(id,:); C=C(:);

id2 = find(log10(C)>-7);
X=X(id2);
Y=Y(id2);
C=C(id2);

figure(1);clf
agutwocolumn(0.6)
wysiwyg
h=scatter(log10(X),log10(Y),[],log10(C),'filled','MarkerFaceAlpha',0.2)
colorbar
grid on
caxis([-8 -3])
xlabel('log_{10}N^2','fontsize',16)
ylabel('log_{10}T_z','fontsize',16)



%% find data points where chipod epsilon is very large

[I,J] = find( log10(chipod.eps)>-4);
unique(chipod.cnum(J));

figure(1);clf
plot(chipod.cnum(J),chipod.P(I),'k.')
axis ij
grid on
xlabel('cast #')
ylabel('P')
title('locations where log_{10}[ \epsilon_{\chi}] > -5')

[660:730]
[2282 2283 2391 2762 2953]

%%

ic=find(chipod.cnum==2762)
figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
plot(log10(chipod.N2(:,ic)),chipod.P)
axis ij
grid on
xlabel('N^2')

subplot(222)
plot(log10(chipod.Tz(:,ic)),chipod.P)
axis ij
grid on
xlabel('T_z')

subplot(223)
plot(log10(chipod.chi(:,ic)),chipod.P)
axis ij
grid on
xlabel('\chi')

subplot(224)
plot(log10(chipod.eps(:,ic)),chipod.P)
axis ij
grid on
xlabel('\epsilon')

%%
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/chipod_method_bin/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_2762_avg.mat')

figure(2);clf
plot(avg.T,avg.P)
axis ij
shg

%%