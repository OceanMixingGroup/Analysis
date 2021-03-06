%~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_tiwe_avg_profiles.m
%
% Combine profiles made in Run_tiwe_AP.m into one structure with common
% depth grid.
%
%----------------
% 2/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

tiwe_patches_paths

Flist = dir( fullfile(path_cham_avg,'*mat') )
Nprof=length(Flist)

% Make empty structure for combined data
cham = struct();
pvec      = 1:1:210;
EmpMat    = nan*ones(length(pvec),Nprof);
cham.T    = EmpMat;
cham.S    = EmpMat;
cham.CHI  = EmpMat;
cham.DTDZ = EmpMat;
cham.N2   = EmpMat;
cham.EPSILON = EmpMat;
cham.cnum = nan*ones(1,Nprof);
cham.yday = nan*ones(1,Nprof);
cham.P    = pvec(:);

hb=waitbar(0,'combining TIWE avg profiles')

for ic=1:length(Flist)
    
    waitbar(ic/Nprof,hb)
    
    clear fname cnum
    
    fname=Flist(ic).name;
    
    cnum = str2num(fname(6:9));       
    
    load( fullfile( path_cham_avg, fname) )
    
    try
        cham.cnum(ic)     = cnum ;
        cham.yday(ic)     = str2num(head.starttime(end-5:end));
        cham.T(:,ic)      = interp1(avg.P,avg.T,pvec) ;
        cham.S(:,ic)      = interp1(avg.P,avg.S,pvec) ;
        cham.N2(:,ic)     = interp1(avg.P,avg.N2,pvec) ;
        cham.DTDZ(:,ic)   = interp1(avg.P,avg.DTDZ,pvec) ;
        cham.CHI(:,ic)    = interp1(avg.P,avg.CHI,pvec) ;
        cham.EPSILON(:,ic)= interp1(avg.P,avg.EPSILON,pvec) ;
    catch
        disp(['Error on profile ' num2str(cnum)])
    end
    
end

delete(hb)

%%

cham.MakeInfo=['Made' datestr(now) ' w/ Combine_tiwe_avg_profiles.m']
save( fullfile('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/','tiwe_1mavg_combined.mat'),'cham')



%%