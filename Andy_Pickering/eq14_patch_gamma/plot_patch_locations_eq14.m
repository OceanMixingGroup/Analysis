function h=plot_patch_locations_eq14(patch_size_min,usetemp,...
    merge_patches,min_sep)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_patch_locations_eq14.m
%
%
%-----------
% 3/17/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%clear ; close all

saveplots=1

eq14_patches_paths

datdir = fullfile(analysis_dir,project_long,'data'); 

% load binned Chameleon data (to plot epsilon etc.)
clear cham
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

h=figure;clf
agutwocolumn(0.75)
wysiwyg

ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon'
ylabel('P')
xlabel('cast #')


% load my patches
clear patches
patches = load_eq14_patches_comb(patch_size_min,usetemp,merge_patches,min_sep)

hold on
plot(patches.cnum,patches.p1,'.','color',0.2*[1 1 1])
%plot(patches.cnum,patches.p2,'.')
xlim([500 2500])
%xlim([3000 3400])
ylim([0 220])

if merge_patches==1
    title(['min OT ' num2str(patch_size_min) 'm, merged'])
else
    title(['min OT ' num2str(patch_size_min) 'm'])
end

%%

if saveplots==1
    
    if merge_patches==1
        fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patch_locs_merged_minsep_' num2str(min_sep*100)]
    else
        fname=[project_short '_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patch_locs' ]
    end
    
    print(fullfile(fig_dir,fname),'-dpng')
end


clear cham patches

%%
