function chb = load_chipod_avg(path_chipod_bin,project_short,Params,cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Loads profile of chipod-method applied to Chameleon data.
%
% INPUT
% path_chipod_bin
% project_short
% Params
% cnum 
%
% OUTPUT
% chb
%
%---------------
% 5/9/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear avg
if strcmp(project_short,'eq14')
    load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[upper(project_short) '_' sprintf('%04d',cnum) '_avg.mat']))
else
    load( fullfile( path_chipod_bin, ['zsm' num2str(Params.z_smooth) 'm_fmax' num2str(Params.fmax) 'Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128'],[project_short '_' sprintf('%04d',cnum) '_avg.mat']))
end
chb = avg;clear avg

%%