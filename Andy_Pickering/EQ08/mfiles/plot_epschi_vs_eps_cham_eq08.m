%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% plot_epschi_vs_eps_cham_eq08.m
%
% Plot eps_chi vs eps for Chameleon data.
%
%
%~~~~~~~~~~~~~~~~~
% 6/13/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; %close all

Params.gamma    = 0.2;
Params.fmax     = 10  ;
Params.z_smooth = 10 ;
Params.resp_corr= 0  ;
Params.fc       = 99 ;

dz = 2 ;
cnums_to_get = 200:2700;

screen_chi= 1 ;
screen_ml = 1 ;
Pmin      = 20 ;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/mixingsoftware/CTD_Chipod/mfiles/
eq08_patches_paths

save_name = [project_short '_screen_chi_' num2str(screen_chi) '_screen_ml_' num2str(screen_ml) '_Pmin_' num2str(Pmin) '_dz_' num2str(dz) '_'  MakeChiPathStr(Params) '.mat']
if exist(fullfile(analysis_dir,project_short,'data',save_name),'file')==2
    load(fullfile(analysis_dir,project_short,'data',save_name))
else
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end

eps_chi = cham.chi .* cham.N2 / 2 / 0.2 ./ (cham.Tz.^2);

figure%(1);clf
agutwocolumn(0.6)
wysiwyg

subplot(1,3,[1 2])
histogram2( real(log10(cham.eps)), real(log10(eps_chi)), 'DisplayStyle','tile','EdgeColor','none','XBinEdges',[-8.5:0.2:-5],'YBinEdges',[-11:0.2:-5])
xlim([-8.5 -5])
ylim([-11 -5])
xvec = linspace(-8.5,-5,100);
hold on
plot(xvec,xvec,'k--')
plot(xvec,xvec+1,'r--')
plot(xvec,xvec-1,'r--')
xlabel('log_{10}[\epsilon]','fontsize',16)
ylabel('log_{10}[\epsilon_{\chi}]','fontsize',16)


subplot(1,3,3)
histogram( real(log10( eps_chi(:) ./ cham.eps(:))),'edgecolor','none','Normalization','pdf')
grid on
freqline( nanmean(real(log10( eps_chi(:) ./ cham.eps(:)))))
xlim([-3.5 2])
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]','fontsize',16)

eq08_patches_paths
figname = [project_short '_epschi_vs_eps']
print(fullfile(fig_dir, figname), '-dpng')



%%