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

clear ; close all

Params.gamma    = 0.2;
Params.fmax     = 10  ;
Params.z_smooth = 1 ;
Params.resp_corr= 0  ;
Params.fc       = 99 ;

dz = 2 ;
cnums_to_get = 200:2700;

screen_chi= 1 ;
screen_ml = 0 ;
Pmin      = 0 ;

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

figure(1);clf
agutwocolumn(0.6)
wysiwyg

histogram2( real(log10(cham.eps)), real(log10(eps_chi)), 'DisplayStyle','tile','EdgeColor','none')
xlim([-8.5 -5])
ylim([-11 -5])
xvec = linspace(-8.5,-5,100);
hold on
plot(xvec,xvec,'k--')
plot(xvec,xvec+1,'r--')
plot(xvec,xvec-1,'r--')
xlabel('\epsilon','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)

eq08_patches_paths
figname = [project_short '_epschi_vs_eps']
print(fullfile(fig_dir, figname), '-dpng')

%%