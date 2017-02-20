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

data_dir = '/Users/Andy/Cruises_Research/ChiPod/TIWE/data/avg/' ;

Flist = dir( fullfile(data_dir,'*mat') )
Nprof=length(Flist)

% Make empty structure for combined data
cham=struct();
pvec=1:1:210;
EmpMat=nan*ones(length(pvec),Nprof);
cham.T=EmpMat;
cham.S=EmpMat;
cham.CHI=EmpMat;
cham.DTDZ=EmpMat;
cham.N2=EmpMat;
cham.EPSILON=EmpMat;
cham.cnum=nan*ones(1,Nprof);
cham.P=pvec(:);

hb=waitbar(0)

for ic=1:length(Flist)
   
    waitbar(ic/Nprof,hb,'working')
    
    clear fname cnum
    
    fname=Flist(ic).name;
   
    cnum = str2num(fname(5:8));
    
    load( fullfile( data_dir, fname) )
    
    try 
    cham.cnum(ic)     = cnum ;
    cham.T(:,ic)      = interp1(avg.P,avg.T,pvec) ;
    cham.S(:,ic)      = interp1(avg.P,avg.S,pvec) ;
    cham.N2(:,ic)     = interp1(avg.P,avg.N2,pvec) ;
    cham.DTDZ(:,ic)   = interp1(avg.P,avg.DTDZ,pvec) ;
    cham.CHI(:,ic)    = interp1(avg.P,avg.CHI,pvec) ;
    cham.EPSILON(:,ic)= interp1(avg.P,avg.EPSILON,pvec) ;
    end
    
end

delete(hb)

%%

cham.MakeInfo=['Made' datestr(now) ' w/ Combine_tiwe_avg_profiles.m']
save( fullfile('/Users/Andy/Cruises_Research/ChiPod/TIWE/data/','tiwe_1mavg_combined.mat'),'cham')


%%

clear ; close all

load( fullfile('/Users/Andy/Cruises_Research/ChiPod/TIWE/data/','tiwe_1mavg_combined.mat') )

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(211);
ezpc(cham.cnum,cham.P,cham.T)
xlabel('cnum')
ylabel('P')
colorbar
title('T')

ax2=subplot(212);
ezpc(cham.cnum,cham.P,cham.S)
xlabel('cnum')
ylabel('P')
colorbar
title('S')

%%

figure(2);clf
histogram2( log10(cham.CHI(:)), log10(cham.EPSILON(:)),'DisplayStyle','tile')
xlabel('log_{10}[\chi]','fontsize',16)
ylabel('log_{10}[\epsilon]','fontsize',16)
colorbar
caxis([0 1000])
title('tiwe 1m avg data')


%%

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(411);
ezpc(cham.cnum,cham.P,real(log10(cham.N2)))
xlabel('cnum')
ylabel('P')
colorbar
caxis([-6 -1])
title('log_{10}[N^2]')

ax2=subplot(412);
ezpc(cham.cnum,cham.P,real(log10(cham.DTDZ)))
xlabel('cnum')
ylabel('P')
colorbar
caxis([-4 0])
title('log_{10}[dT/dz]')

ax3=subplot(413);
ezpc(cham.cnum,cham.P,log10(cham.CHI))
xlabel('cnum')
ylabel('P')
colorbar
caxis([-11 -4])
title('log_{10}[\chi]')

ax4=subplot(414);
ezpc(cham.cnum,cham.P,log10(cham.EPSILON))
xlabel('cnum')
ylabel('P')
colorbar
caxis([-11 -5])
title('log_{10}[\epsilon]')

linkaxes([ax1 ax2 ax3 ax4])

%%

fig_dir='/Users/Andy/Cruises_Research/ChiPod/TIWE/figures'
print(fullfile(fig_dir,'tiwe_avgCombine_N2_dtdz_chi_eps'),'-dpng')

SetNotesFigDir
print(fullfile(NotesFigDir,'tiwe_avgCombine_N2_dtdz_chi_eps'),'-dpng')

%% Compute gamma from these values

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

gam=ComputeGamma(cham.N2(:),cham.DTDZ(:),cham.CHI(:),cham.EPSILON(:));

figure(2);clf
histogram(real(log10(gam(:))),'EdgeColor','none','Normalization','pdf')
freqline(log10(0.2))
xlim([-4 2])
grid on
xlabel('log_{10}[\Gamma]')
ylabel('pdf')
title('tiwe 1m bin processed by AP')

%%
fig_dir='/Users/Andy/Cruises_Research/ChiPod/TIWE/figures'
print(fullfile(fig_dir,'tiwe_avgCombine_gamma'),'-dpng')

SetNotesFigDir
print(fullfile(NotesFigDir,'tiwe_avgCombine_gamma'),'-dpng')

%%