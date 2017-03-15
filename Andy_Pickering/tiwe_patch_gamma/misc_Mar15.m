
%%

patches = load_tiwe_patches_comb(0.4,1,0.0)
id = find(patches.R2>0.5) ;

figure(2);clf
histogram(log10(patches.gam_line),'Normalization','pdf')
hold on
histogram(log10(patches.gam_line(id)),'Normalization','pdf')

%%

figure(1);clf
plot(patches.R2,patches.gam_line,'.')

%%