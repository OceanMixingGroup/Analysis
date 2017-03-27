%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% compare_patches_diff_experiments
%
% Compare Patch Distributions from tiwe,eq08,eq14?
%
%
% 3/23/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.4;
usetemp = 1;
merge_patches = 0;
min_sep = 0;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

patches_tiwe = load_patches_comb('tiwe',patch_size_min, usetemp, merge_patches, min_sep)
patches_eq08 = load_patches_comb('eq08',patch_size_min, usetemp, merge_patches, min_sep)
patches_eq14 = load_patches_comb('eq14',patch_size_min, usetemp, merge_patches, min_sep)


depth_range = [80 200]
id_tiwe = find(patches_tiwe.p1>depth_range(1) & patches_tiwe.p2<depth_range(2));
id_eq08 = find(patches_eq08.p1>depth_range(1) & patches_eq08.p2<depth_range(2));
id_eq14 = find(patches_eq14.p1>depth_range(1) & patches_eq14.p2<depth_range(2));


%%

figure(4);clf
agutwocolumn(1)
wysiwyg

xl=[-3 1.5];

ax1=subplot(311)
h1 = histogram(real(log10(patches_tiwe.gam_line(id_tiwe))),'EdgeColor','none','Normalization','pdf');
freqline(log10(0.2))
xlim(xl)
grid on
ylabel('pdf')
SubplotLetterMW(['tiwe, N=' num2str(length(id_tiwe))])

ax2=subplot(312)
h2 = histogram(real(log10(patches_eq08.gam_line(id_eq08))),h1.BinEdges,'EdgeColor','none','Normalization','pdf');
xlim(xl)
freqline(log10(0.2))
%freqline(log10(0.1),'r--')
%freqline(log10(0.3),'r--')
grid on
ylabel('pdf')
SubplotLetterMW(['eq08, N=' num2str(length(id_eq08))])


ax3=subplot(313)
h3 = histogram(real(log10(patches_eq14.gam_line(id_eq14))),h1.BinEdges,'EdgeColor','none','Normalization','pdf');
xlim(xl)
freqline(log10(0.2))
grid on
ylabel('pdf')
%SubplotLetterMW('eq14')
SubplotLetterMW(['eq14, N=' num2str(length(id_eq14))])


linkaxes([ax1 ax2 ax3],'x')
xlabel('log_{10}[\gamma line]','fontsize',16)

%%

fig_dir = '/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/figures/'
ChkMkDir(fig_dir)
figname =  ['tiwe_eq08_eq14_gam_compare_minOT_' num2str(100*patch_size_min)...
    '_usetemp_' num2str(usetemp) '_zrange_' num2str(depth_range(1)) '_'...
    num2str(depth_range(2))]
print( fullfile( fig_dir,figname), '-dpng')

%%

figure(3);clf
agutwocolumn(1)
wysiwyg

xl=[-3 1.5];

ax1=subplot(311)
h1 = histogram(real(log10(patches_tiwe.gam_bulk(id_tiwe))),'EdgeColor','none','Normalization','pdf');
freqline(log10(0.2))
xlim(xl)
grid on
ylabel('pdf')
SubplotLetterMW('tiwe')

ax2=subplot(312)
h2 = histogram(real(log10(patches_eq08.gam_bulk(id_eq08))),h1.BinEdges,'EdgeColor','none','Normalization','pdf');
xlim(xl)
freqline(log10(0.2))
%freqline(log10(0.1),'r--')
%freqline(log10(0.3),'r--')
grid on
ylabel('pdf')
SubplotLetterMW('eq08')

ax3=subplot(313)
h3 = histogram(real(log10(patches_eq14.gam_bulk(id_eq14))),h1.BinEdges,'EdgeColor','none','Normalization','pdf');
xlim(xl)
freqline(log10(0.2))
grid on
ylabel('pdf')
SubplotLetterMW('eq14')

linkaxes([ax1 ax2 ax3],'x')
xlabel('log_{10}[\gamma bulk]','fontsize',16)

%% Plot patch locations over chameleon epsilon for all 3 experiments

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(311)



%%