
%%

clear ; close all

Params.gamma = 0.2 ;
Params.fmax  = 7 ;
Params.z_smooth = 10 ;

%dz = 10 % bin size
zmin=0  ;
zmax=200;

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/

eq14_patches_paths
Pmin = 80;
%cnums_to_get = get_cham_cnums_eq14 ;
cnums_to_get = 2000:3000;


figure(2);clf
agutwocolumn(1)
wysiwyg

iax=1

for dz=[1 10 30]
    
    clear chipod cham
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,zmin,zmax)
    
    % Nan values in mixed layer
    clear ib
    ib = find(cham.P<Pmin);
    cham.eps(ib) = nan;
    
    clear ib
    ib = find(chipod.P<Pmin);
    chipod.eps(ib) = nan;
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(cham.eps./cham.N2)),log10(cham.chi./(cham.Tz.^2)),80,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(1e-7,1e-1,100);
    h1=plot( log10(xvec), log10(xvec*2*0.2),'k-');
    h2=plot( log10(xvec), log10(xvec*2*0.1),'r-');
    h3=plot( log10(xvec), log10(xvec*2*0.05),'c-');
    ylim([-7.5 -1])
    xlim([-5.5 -1])
    ylabel('log_{10} [\chi / T_{z}^{2}]','fontsize',16)
    if iax>4
        xlabel('log_{10} [\epsilon / N^{2}]','fontsize',16)
    end
    legend([h1 h2 h3],['\gamma=0.2'],['\gamma=0.1'],['\gamma=0.05'],'location','best')
        title([project_short ' Chameleon ' num2str(dz) 'm binned, >' num2str(Pmin) 'db'])
    
    iax=iax+1
    
    subplot(3,2,iax)
    hh=histogram2(  real(log10(cham.eps)),log10(chipod.eps),80,'DisplayStyle','tile')
    grid on
    hold on
    xvec=linspace(-11,-4,100);
    plot(xvec,xvec,'k--')
    ylim([-11 -4])
    xlim([-11 -4])
    ylabel('log_{10} [\epsilon_{\chi}]','fontsize',16)
    if iax>4
        xlabel('log_{10} [\epsilon ]','fontsize',16)
    end
    
    iax=iax+1
    
end

%%
figname=['eq14_NormScat_chiVscham_diff_dz']
print(fullfile(fig_dir,figname),'-dpng')

%%