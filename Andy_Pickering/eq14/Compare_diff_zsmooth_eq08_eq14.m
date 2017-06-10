%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_diff_zsmooth_eq08_eq14.m
%
% Compare chipod estimates, eps_chi vs eps for different values of z_smooth
%
%
% - ComputeChi_Chameleon_Eq08.m
% - 
%
%----------------------
% 5/17/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

%project = 'eq08'
project = 'eq14'

%cnums_to_get = 200:2800%
cnums_to_get = get_cham_cnums_eq14 ;

if strcmp(project','eq14')
bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
cnums_to_get = setdiff(cnums_to_get,bad_prof);
end

eval([project '_patches_paths'])
%eq14_patches_paths

screen_chi = 1 ;
Pmin       = 20;
screen_ml  = 1 ;

Params.gamma    = 0.2;
Params.fmax     = 7 ;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/code

dz = 2 ;

Params.z_smooth = 1 ;
[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

% Get data again for different fmax value
Params.z_smooth = 10 ;
[chipod2, cham2] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);


%% histograms of chipod chi,epsilon 

figure(1);clf

subplot(211)
h1 = histogram(log10(chipod.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
hold on
h2 = histogram(log10(chipod2.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
% h3 = histogram(log10(chipod3.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
% if strcmp(project,'eq14')
%     h4 = histogram(log10(chipod4.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
%     legend([h1 h2 h3 h4],'32hz','15hz','7hz','10hz')
% else
%     legend([h1 h2 h3],'32hz','15hz','10hz')
% end
     legend([h1 h2 ],'1m','10m')
     
grid on
xlim([-11 -3])
xlabel('log_{10}[\chi_{\chi}]')
ylabel('pdf')

subplot(212)
h1 = histogram(log10(chipod.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
hold on
h2 = histogram(log10(chipod2.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
% h3 = histogram(log10(chipod3.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
% if strcmp(project,'eq14')
%     h4 = histogram(log10(chipod4.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
%     legend([h1 h2 h3 h4],'32hz','15hz','7hz','10hz')
% else
%    legend([h1 h2 h3],'32hz','15hz','7hz')
%end

     legend([h1 h2 ],'1m','10m')


grid on
xlim([-8.5 -4])
xlabel('log_{10}[\epsilon_{\chi}]')

%
figname = [project '_chi_eps_histograms_diff_zsmooth']
%print( fullfile(fig_dir, figname), '-dpng')

%% histograms of ratio of chipod:cham chi,epsilon 

figure(1);clf

subplot(211)
h1 = histogram(log10(chipod.chi ./ cham.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
hold on
h2 = histogram(log10(chipod2.chi ./ cham.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
% h3 = histogram(log10(chipod3.chi ./ cham.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
% if strcmp(project,'eq14')
%     h4 = histogram(log10(chipod4.chi ./ cham.chi),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
%     legend([h1 h2 h3 h4],'32hz','15hz','7hz','10hz')
% else
%     legend([h1 h2 h3],'32hz','15hz','10hz')
% end

grid on
xlim([-3 3])
xlabel('log_{10}[\chi_{\chi}/\chi]')
ylabel('pdf')
     legend([h1 h2 ],'1m','10m')
     
subplot(212)
h1 = histogram(log10(chipod.eps ./ cham.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
hold on
h2 = histogram(log10(chipod2.eps ./ cham.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
% h3 = histogram(log10(chipod3.eps ./ cham.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2)
% if strcmp(project,'eq14')
%     h4 = histogram(log10(chipod4.chi ./ cham.eps),'Normalization','pdf','DisplayStyle','stair','LineWidth',2);
%     legend([h1 h2 h3 h4],'32hz','15hz','7hz','10hz')
% else
%     legend([h1 h2 h3],'32hz','15hz','7hz')
% end

grid on
xlim([-3 3])
xlabel('log_{10}[\epsilon_{\chi}/\epsilon]')
     legend([h1 h2 ],'1m','10m')
     
figname = [project '_chi_eps_ratio_histograms_diff_fmax']
%print( fullfile(fig_dir, figname), '-dpng')

%% Plot 2D histograms chipod vs cham values of chi and epsilon

figure(3);clf
agutwocolumn(1)
wysiwyg

if strcmp(project,'eq14')
    rr=2;cc=2;
else
rr=2;cc=2;
end

subplot(rr,cc,1)
histogram2( log10(cham.chi(:)), log10(chipod.chi(:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-12 -4])
ylim([-12 -4])
xlabel('\chi','fontsize',16)
ylabel('\chi_{\chi}','fontsize',16)
title('1m')

subplot(rr,cc,2)
histogram2( log10(cham.eps(:)), log10(chipod.eps(:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-8.5 -4])
ylim([-8.5 -4])
xlabel('\epsilon ','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)
title('1m')
caxis([0 550])

subplot(rr,cc,3)
histogram2( log10(cham.chi(:)), log10(chipod2.chi(:)), 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-12 -4])
ylim([-12 -4])
xlabel('\chi','fontsize',16)
ylabel('\chi_{\chi}','fontsize',16)
title('10m')

subplot(rr,cc,4)
histogram2( log10(cham.eps(:)), log10(chipod2.eps(:)),50, 'DisplayStyle','tile')
hold on
xvec=linspace(-11,-4,100);
plot(xvec,xvec,'k--')
plot(xvec,xvec-1,'r--')
plot(xvec,xvec+1,'r--')
xlim([-8.5 -4])
ylim([-8.5 -4])
xlabel('\epsilon ','fontsize',16)
ylabel('\epsilon_{\chi}','fontsize',16)
title('10m')
caxis([0 550])
%


%%
