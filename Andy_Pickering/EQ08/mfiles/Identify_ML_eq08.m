%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Identify_ML_eq08.m
%
% Find convective regions near surface in EQ08, to be excluded in analysis
%
%
%-----------------
% 5/19/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/

eq08_patches_paths
%
cnums_to_get = 200:2700;

Params.gamma    = 0.2;
Params.fmax     = 15  ;
Params.z_smooth = 10;
Params.resp_corr= 0;
Params.fc       = 99;

dz = 1 ;
Pmin=0;
screen_chi=0;
screen_ml=0;

[chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);

%

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(211);
ezpc(cham.cnum,cham.P,real(log10(cham.N2)))
caxis([-6 -2])
cb=colorbar
ylim([0 150])
cb.Label.String='log_{10}N^2';cb.FontSize=16;
hline(80,'--')

ax2=subplot(212);
ezpc(cham.cnum,cham.P,log10(cham.eps))
caxis([-11 -4])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';cb.FontSize=16;
ylim([0 150])
hline(80,'--')

linkaxes([ax1 ax2])

%

thresh=0.04
zml=nan*ones(1,length(1:2700));
for i=1:2700%length(cham.cnum)
    
    try
        
        clear avg
        
        cnum = i
        
        % load profile
        load( fullfile( path_cham_avg,['eq08_' sprintf('%04d',cnum) '_avg.mat']) )
        
        clear ig sigma_surf
        ig = find(~isnan(avg.SIGMA(:)));
        sigma_surf = avg.SIGMA(ig(1));
        id=find( avg.SIGMA(:)>(sigma_surf+thresh) );
        if ~isempty(id)
            zml(i)=avg.P(id(1));
        end
    end
end

hold on

zml_cnum=1:2700;

plot(zml_cnum,zml,'k')
subplot(211)
hold on
plot(zml_cnum,zml,'k')

linkaxes([ax1 ax2])

%% save file with depth for each cast to use in thresholding

save( fullfile(analysis_dir,project_short,'data','eq08_mldepths.mat'),'zml','zml_cnum')

%%