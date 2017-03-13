
%% compare # patches w/ different params (all yday 324-327)

clear ; close all

patch_size_min = 1  % min patch size
usetemp = 1
merge_patches = 0 ;
min_sep = 0.15 ;

patches=load_tiwe_patches_comb(patch_size_min, usetemp, merge_patches, min_sep)

%
%id=isin(patches.yday,[323 327])
id = find(patches.p1>60 & patches.p2<200);
length(id)

%% Plot hist of gamma for different patch params

p1=load_tiwe_patches_comb(0.4, 1, 0, .15)
id1 = find(p1.p1>60 & p1.p2<200);

p2=load_tiwe_patches_comb(0.4, 0, 0, .15)
id2 = find(p2.p1>60 & p2.p2<200);

p3=load_tiwe_patches_comb(1, 1, 0, .15)
id3 = find(p3.p1>60 & p3.p2<200);

p4=load_tiwe_patches_comb(1, 0, 0, .15)
id4 = find(p4.p1>60 & p4.p2<200);

figure(1);clf
histogram( real(log10(p1.gam_bin(id1))),'Normalization','pdf','Edgecolor','none')
hold on
histogram( real(log10(p2.gam_bin(id2))),'Normalization','pdf','Edgecolor','none')
histogram( real(log10(p3.gam_bin(id3))),'Normalization','pdf','Edgecolor','none')
histogram( real(log10(p4.gam_bin(id4))),'Normalization','pdf','Edgecolor','none')
xlim([-4 2])
freqline(log10(0.2))
grid on

figure(2);clf
h1 = histogram( real(log10(p1.gam_line(id1))),'Normalization','pdf','Edgecolor','none')
hold on
h2 = histogram( real(log10(p2.gam_line(id2))),'Normalization','pdf','Edgecolor','none')
h3 = histogram( real(log10(p3.gam_line(id3))),'Normalization','pdf','Edgecolor','none')
h4 = histogram( real(log10(p4.gam_line(id4))),'Normalization','pdf','Edgecolor','none')
xlim([-4 2])
freqline(log10(0.2))
grid on
legend([h1 h2 h3 h4],'40cm, temp','40cm, dens','1m temp','1mdens')

