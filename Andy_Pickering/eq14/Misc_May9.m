%
%
% Compare Get_binned_data_avg_profile_v2 to Get_all_chipod_cham_data w/
% averaging done after
%
%
%
%%

clear ; close all

Params.gamma = 0.2;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

dz=10; % bin size

eq14_patches_paths

plot_dir = fullfile(fig_dir,'chi_eps_profiles_wPoints')
ChkMkDir(plot_dir)

Pmin=20;
screen_chi=1
screen_ml=0

eps_cham = [];
eps_chi = [];
eps_bias = [];
P_all = [];
N2_all = [];
Tz_all = [];

for cnum=1000%1:1:3000
    
    clear cnum_range
    cnum_range = [cnum cnum+50];
    
    clear cnums_to_get
    cnums_to_get = [cnum_range(1) : cnum_range(2) ];
    bad_prof=[2282 2283 2391 2762 2953]; % profiles where temp. is bad
    cnums_to_get = setdiff(cnums_to_get,bad_prof);
    
    clear chipod cham
    [chipod, cham] = Get_all_chipod_cham_data(path_chipod_bin,...
        path_cham_avg,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml);
    %
    
    ib=find(log10(chipod.eps)>-5);
    chipod.eps(ib)=nan;
    ib=find(log10(cham.eps)>-5);
    % cham.eps(ib)=nan;
    
    [chi_cham_avg, p_cham, ~] = binprofile(cham.chi, cham.P, 0, dz, 200,1);
    [chi_chi_avg, p_chi, ~] = binprofile(chipod.chi, chipod.P, 0, dz, 200,1);
    
    [eps_cham_avg, p_cham, ~] = binprofile(cham.eps, cham.P, 0, dz, 200,1);
    [eps_chi_avg, p_chi, ~] = binprofile(chipod.eps, chipod.P, 0, dz, 200,1);
    
    %~~~ do same with Get_binned_data_avg_profile_v2.m
    clear chipod cham
    [chipod, cham] = Get_binned_data_avg_profile_v2(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,Pmin,screen_chi,screen_ml)
    
end

% plot: should be exactly the same?

figure(1);clf
loglog(cham.eps,eps_cham_avg,'o')
hold on
loglog(chipod.eps,eps_chi_avg,'d')
grid on
xvec=linspace(1e-9,1e-4,100);
plot(xvec,xvec,'k--')

%%