%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_ChipodMethodPatches.m
%
% This script combines (1) chi-pod estimates at patches using constant
% gamma, (2) chi-pod estimates at patches using actual gammas , (3) chi-pod
% estimates using binned data, and (4) chameleon estimates at patches, for
% all profiles and patches so we can compare them etc.
%
% - chi-pod estimates for patches are made in : ComputeChi_Chameleon_Eq14_PATCHES.m
%   - Chameleon chi/eps values at patch locations are also saved as
%   'eps_bin', 'eps_patch' etc
% - chi-pod estimates for binned profiles are made w/ ComputeChi_Chameleon_Eq14.m
%
% Produces a structure 'AllEps'
%
% Used to be part of CompareProfiles_bin_patch_cham.m
%
%--------------
% 1/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Combine results from all profiles (so we can make scatter plot etc) and save

clear ; close all

whN2dTdz = 'line'
%whN2dTdz = 'line_fit'
%whN2dTdz = 'bulk'
Params.gamma = 0.2;
Params.fmax=7

% patch parameters
patch_size_min = 0.4
usetemp = 1

minR2 = 0.5

savedata = 1;

% pre-allocate empty arrays
P=[];

eps_patchN2dTdz_constGam = [] ;
eps_patchN2dTdzGam    = [] ;
eps_cham_bin  = [] ;
eps_cham_patch= [] ;
eps_chipod_binned    = [] ;

chi_patchN2dTdz_constGam = [] ;
chi_patchN2dTdzGam    = [] ;
chi_cham_bin  = [] ;
chi_cham_patch= [] ;
chi_chipod_binned    = [] ;

eq14_patches_paths

dir1 = fullfile(analysis_dir,project_long,'data','ChipodPatches')

hb=waitbar(0)
%
for cnum=1:3100
    waitbar(cnum/3100,hb)
    try
        
        % patch N^2,dTdz w/ constant gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdz_constGam=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_minR2_' num2str(minR2)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdzGam=avg;clear avg
        
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_chipod_binned=avg;clear avg
        
        % get these (binned) values at patch locations
        clear eps3 chi3 pp
        eps3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        chi3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        
        for ip=1:length(avg_patchN2dTdz_constGam.P)
            [val,I]=nanmin(abs(avg_chipod_binned.P-avg_patchN2dTdz_constGam.P(ip)));
            eps3(ip)=avg_chipod_binned.eps1(I);
            chi3(ip)=avg_chipod_binned.chi1(I);
        end
        
        P = [P(:) ; avg_patchN2dTdz_constGam.P(:) ];
        
        % regular chi-pod method epsilons (1m smooth, gam=0.2 etc)
        eps_chipod_binned  = [eps_chipod_binned   ; eps3(:)        ];
        chi_chipod_binned  = [chi_chipod_binned   ; chi3(:)        ];
        
        % chi-pod estimates at patches, using constant gamma
        eps_patchN2dTdz_constGam = [eps_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.eps1(:) ];
        chi_patchN2dTdz_constGam = [chi_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.chi1(:) ];
        
        % chi-pod estimates at patches, using actual patch gammas
        eps_patchN2dTdzGam       = [eps_patchN2dTdzGam       ; avg_patchN2dTdzGam.eps1(:)       ];
        chi_patchN2dTdzGam       = [chi_patchN2dTdzGam       ; avg_patchN2dTdzGam.chi1(:)       ];
        
        % chameleon epsilon binned data at patch locations
        eps_cham_bin = [eps_cham_bin ; avg_patchN2dTdz_constGam.eps_bin(:)];
        chi_cham_bin = [chi_cham_bin ; avg_patchN2dTdz_constGam.chi_bin(:)];
        
        % chameleon epsilon values computed over patches only
        eps_cham_patch = [eps_cham_patch ; avg_patchN2dTdzGam.eps_patch(:) ] ;
        chi_cham_patch = [chi_cham_patch ; avg_patchN2dTdzGam.chi_patch(:) ] ;
        
    end % try
    
end % cnum
delete(hb)
%
ib = find(log10(eps_cham_bin)<-8.5);
eps_cham_bin(ib) = nan;

AllEps = struct('eps_chipod_binned',eps_chipod_binned,'eps_patchN2dTdz_constGam',...
    eps_patchN2dTdz_constGam,'eps_patchN2dTdzGam',eps_patchN2dTdzGam,...
    'eps_cham_bin',eps_cham_bin,'eps_cham_patch',eps_cham_patch,'chi_chipod_binned',chi_chipod_binned,'chi_patchN2dTdz_constGam',...
    chi_patchN2dTdz_constGam,'chi_patchN2dTdzGam',chi_patchN2dTdzGam,...
    'chi_cham_bin',chi_cham_bin,'chi_cham_patch',chi_cham_patch,'P',P)
AllEps.MakeInfo = ['Made ' datestr(now) ' w/ Combine_ChipodMethodPatches.m']

if savedata==1
    % save data
    sav_name = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_minOT' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_chipodmethods.mat'] ;
    sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches'
    save(fullfile(sav_dir,sav_name),'AllEps')
end

%%


%id = find(AllEps.P>0 & AllEps.P<200);
id = find(AllEps.P>60 & AllEps.P<200);

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
histogram2( log10(AllEps.chi_cham_bin(:)), log10(AllEps.chi_cham_patch(:)), 50, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\chi cham bin','fontsize',16)
ylabel('\chi cham patch','fontsize',16)

subplot(222)
histogram2( log10(AllEps.chi_cham_patch(id)), log10(AllEps.chi_patchN2dTdzGam(id)),40, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\chi cham patch','fontsize',16)
ylabel('\chi chipod actual gam','fontsize',16)

subplot(223)
histogram2( log10(AllEps.chi_cham_patch(id)), log10(AllEps.chi_chipod_binned(id)),40, 'DisplayStyle','Tile')
xlim([-10 -4])
ylim([-10 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\chi cham patch','fontsize',16)
ylabel('\chi chipod binned','fontsize',16)

subplot(224)
histogram2( log10(AllEps.chi_cham_patch(id)), log10(AllEps.chi_patchN2dTdz_constGam(id)),40, 'DisplayStyle','Tile')
xlim([-10 -4])
ylim([-10 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\chi cham patch','fontsize',16)
ylabel('\chi chipod patch','fontsize',16)

%
eq14_patches_paths
figname = ['ChiPatchChiCompare_N2dTdz_' num2str(whN2dTdz) '_minOT' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) ] ;
print(fullfile(fig_dir,figname),'-dpng')


%%

%id = find(AllEps.P>0 & AllEps.P<200);
id = find(AllEps.P>60 & AllEps.P<200);

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
histogram2( log10(AllEps.eps_cham_bin(:)), log10(AllEps.eps_cham_patch(:)), 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\epsilon cham bin','fontsize',16)
ylabel('\epsilon cham patch','fontsize',16)

subplot(222)
histogram2( log10(AllEps.eps_cham_patch(id)), log10(AllEps.eps_patchN2dTdzGam(id)),50, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\epsilon cham patch','fontsize',16)
ylabel('\epsilon chipod actual gam','fontsize',16)

subplot(223)
histogram2( log10(AllEps.eps_cham_patch(id)), log10(AllEps.eps_chipod_binned(id)),50, 'DisplayStyle','Tile')
xlim([-10 -4])
ylim([-10 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\epsilon cham patch','fontsize',16)
ylabel('\epsilon chipod binned','fontsize',16)

subplot(224)
histogram2( log10(AllEps.eps_cham_patch(id)), log10(AllEps.eps_patchN2dTdz_constGam(id)),50, 'DisplayStyle','Tile')
xlim([-10 -4])
ylim([-10 -4])
hold on
xvec=linspace(-12,-3,100);
loglog(xvec,xvec,'k--')
xlabel('\epsilon cham patch','fontsize',16)
ylabel('\epsilon chipod patch','fontsize',16)

%%

eq14_patches_paths
figname = ['ChiPatchEpsCompare_N2dTdz_' num2str(whN2dTdz) '_minOT' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) ] ;
print(fullfile(fig_dir,figname),'-dpng')

%% histogram of ratio of different epsilon estimates to actual epsilon

figure(2);clf
agutwocolumn(0.6)
wysiwyg

hbin = histogram(log10(AllEps.eps_chipod_binned(id)./AllEps.eps_cham_patch(id)),'EdgeColor','none','Normalization','pdf');
hold on
hpatch = histogram(log10(AllEps.eps_patchN2dTdz_constGam(id) ./AllEps.eps_cham_patch(id)),'EdgeColor','none','Normalization','pdf');
%histogram(log10(AllEps.eps_patchN2dTdzGam(id) ./AllEps.eps_cham_patch(id)),'EdgeColor','none','Normalization','pdf')
xlim([-5 2])
grid on
freqline(log10(1))
legend([hbin hpatch],'bin','patch')
xlabel('log_{10}[\epsilon_{\chi}/\epsilon_{\epsilon}]','fontsize',16)
ylabel('pdf','fontsize',16)

%
eq14_patches_paths
figname = ['ChiPatchEpsHist_N2dTdz_' num2str(whN2dTdz) '_minOT' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) ]
print(fullfile(fig_dir,figname),'-dpng')

%%