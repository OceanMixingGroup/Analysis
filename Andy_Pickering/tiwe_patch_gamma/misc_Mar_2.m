
%%

figure(1);clf
id=find(patches.Lt<5);
histogram(log10(patches.gam_line(:)),'Normalization','pdf')
hold on
histogram(log10(patches.gam_line(id)),'Normalization','pdf')
freqline(0.2)

%%

figure(2);clf
loglog(patches.gam_line,patches.Lt,'.')
freqline(0.2)
xlim([1e-2 1e1])
xlabel('\gamma')
ylabel('L_t')

%% compare merged/unmerged

clear ; close all

day_range=[307 329]% all profiles
%day_range=[324 327]% ydays in Smyth et al


load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_cham_minOT_15_usetemp_1_patches_diffn2dtdzgamma_merged_minsep_15.mat')
p1=patches; clear patches
id1=find(p1.yday>=day_range(1) & p1.yday<=day_range(2));


load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_cham_minOT_15_usetemp_1_patches_diffn2dtdzgamma.mat')
p2=patches; clear patches
id2=find(p2.yday>=day_range(1) & p2.yday<=day_range(2));


figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(311)
hmerg=histogram(real(log10(p1.gam_bin(id1))),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_bin(id2))),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')
legend([hmerg h2],'merged','not','location','best')
grid on
title('binned')
ylabel('pdf')

%figure(2);clf
subplot(312)
hmerg=histogram(real(log10(p1.gam_line(id1))),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_line(id2))),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')
legend([hmerg h2],'merged','not','location','best')
grid on
title('line')
ylabel('pdf')

%figure(3);clf
subplot(313)
hmerg=histogram(real(log10(p1.gam_bulk(id1))),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_bulk(id2))),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 1])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')
legend([hmerg h2],'merged','not','location','best')
grid on
title('bulk')
ylabel('pdf')

%% compare different min OT sizes

% ** No difference in gamma between minOT=15,25cm?

clear ; close all

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_cham_minOT_15_usetemp_1_patches_diffn2dtdzgamma.mat')
p1=patches; clear patches

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/tiwe_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')
p2=patches; clear patches

figure(1);clf
h1=histogram(real(log10(p1.gam_bin)),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_bin)),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 2])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')

figure(2);clf
h1=histogram(real(log10(p1.gam_bulk)),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_bulk)),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 2])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')

figure(3);clf
h1=histogram(real(log10(p1.gam_line)),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(p2.gam_line)),'Normalization','pdf','EdgeColor','none')
xlim([-3.5 2])
freqline(log10(0.2))
xlabel('log_{10}[\gamma]')

%%