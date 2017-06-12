%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareOldNewEQ08Cham.m
%
% Compare the chameleon data I processed (using smaller famx) to the
% processed EQ08 data I got (processed by Sasha?)
%
% Want to see how it changes chi and eps, and make sure I didn't totally
% screw something up when re-processing.
%
% 6/12/17 -  A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Volumes/SP PHD U3/NonBackup/EQ08/processed/eq08_sum.mat')
cham_orig=cham ; clear cham

Params.gamma    = 0.2;
Params.fmax     = 10  ;
Params.z_smooth = 10 ;
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
if exist(fullfile(analysis_dir,project_short,'Data',save_name),'file')==2
    load(fullfile(analysis_dir,project_short,'Data',save_name))
else
    [chipod, cham] = Get_and_bin_profiles(path_chipod_bin,path_cham_avg,dz,Params,cnums_to_get,project_short,0,200,Pmin,screen_chi,screen_ml);
end


figure(1);clf
agutwocolumn(0.8)
wysiwyg

subplot(211)
h0=histogram(log10(cham_orig.CHI),'Normalization','pdf');
hold on
h1=histogram(log10(cham.chi),'Normalization','pdf');
xlim([-12 -3])
grid on
legend([h0 h1],'orig','reproc')
xlabel('log_{10}[\chi]')
ylabel('pdf')

ib=find(log10(cham_orig.EPSILON)<-8.5);
cham_orig.EPSILON(ib)=nan;

subplot(212)
histogram(log10(cham_orig.EPSILON),'Normalization','pdf')
hold on
histogram(log10(cham.eps),'Normalization','pdf')
grid on
xlabel('log_{10}[\epsilon]')
ylabel('pdf')

%%