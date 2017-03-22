%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% find_good_TS_tiwe_patches.m
%
% Examine T-S relationship for patches and try to find ones with clear T-S
% relationship (do linear fit and use R^2?)
%
%
% 3/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch parameters
patch_size_min=0.4
usetemp=1;
% option to use merged patches
merge_patches = 0 ;
min_sep = 0.15 ;

ot_dir=['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)];

% set paths
tiwe_patches_paths

% load combined patch data for all profiles
patches=load_tiwe_patches_comb(patch_size_min, usetemp, 0, .15) ;
%
addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/


patches.R2 = nan*ones(size(patches.cnum));

for cnum=2836:3711 %1:4000
    
    try
        
        % find patches for this profile
        clear igc Npatches EmpVec
        
        igc = find(patches.cnum==cnum);
        
        Npatches = length(igc);
        % EmpVec=nan*ones(Npatches,1);
        
        % load raw chameleon cast
        clear cal cal2 head
        cal = load_cal_tiwe(cnum);
        cnum_loaded = cnum;
        
        % compute pot. temp, pot. density etc.
        clear s t p lat ptmp sgth
        s = cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
        %s=cal.SAL_sm(1:end-1);
        %s = smooth(s,20);
        t = cal.T1 (1:end-1);
        p = cal.P  (1:end-1);
        
        for ip=1:Npatches
            
            clear x y idz
            idz = isin(cal.P,[patches.p1(igc(ip)) patches.p2(igc(ip))]);
            x=cal.T1(idz);
            y=cal.SAL(idz);
            [P,S] = polyfit(x,y,1);
            
            makeplots = 0 ;
            if makeplots==1
                figure(1);clf
                
                subplot(121)
                plot(cal.T1(idz),cal.P(idz))
                axis ij
                grid on
                xlabel('T')
                xlabel('P')
                
                subplot(122)
                plot(x,y,'o')
                grid on
                hold on
                plot(x,polyval(P,x))
                xlabel('T')
                ylabel('S')
                
            end
            
            % compute R^2 for fit
            R2 = compute_R2(x,y,P) ;
            
            patches.R2(igc(ip))=R2;
            %pause
            
        end % each patch
        
    end % try
    
end % cnum

%%

figure(1);clf
histogram(patches.R2)

%%

N=length(patches.R2)
id=find(patches.R2>0.6);

length(id)/N

%%

figure(1);clf
plot(patches.cnum,patches.p1,'.')
axis ij
hold on
plot(patches.cnum(id),patches.p1(id),'rd')
xlabel('cast #')
ylabel('p[db]')




%%

addpath /Users/Andy/Cruises_Research/Analysis/Andy_Pickering/gen_mfiles/
gam=ComputeGamma(patches.n2_line,patches.dtdz_line,patches.chi,patches.eps);

id=find(patches.p1>60 & patches.p2<200 & patches.R2>0.5);

figure(2);clf
histogram(log10(gam),'Normalization','pdf')
hold on
histogram(log10(gam(id)),'Normalization','pdf')
freqline(log10(0.2))

%%

igc = find(patches.R2>0.75);

patches_good      = struct();
patches_good.p1   = patches.p1(igc) ;
patches_good.p2   = patches.p2(igc) ;
patches_good.chi  = patches.chi(igc) ;
patches_good.eps  = patches.eps(igc) ;
patches_good.cnum = patches.cnum(igc) ;
patches_good.R2   = patches.R2(igc) ;

ip=410

cnum = patches_good.cnum(ip)
p1= patches_good.p1(ip)
p2= patches_good.p2(ip)

% load raw chameleon cast
clear cal cal2 head
cal = load_cal_tiwe(cnum)
cnum_loaded = cnum;

% get data in patch
clear idz
idz = isin(cal.P,[patches_good.p1(ip) patches_good.p2(ip) ]);
t = cal.T1(idz)  ;
p = cal.P(idz)   ;
s = cal.SAL(idz) ;

[P,S] = polyfit(t,s,1);

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
plot(t,p)
axis ij
grid on
xlabel('T')
ylabel('P')

subplot(222)
hraw=plot(s,p);
hold on
hfit = plot(polyval(P,t),p,'linewidth',2);
axis ij
grid on
xlabel('S')
ylabel('P')
legend([hraw hfit],'raw','fit','location','best')

subplot(223)
hraw=plot(t,s,'o');
grid on
hold on
hfit=plot(t,polyval(P,t));
xlabel('T')
ylabel('S')
legend([hraw hfit],'raw','fit','location','best')

sgth=sw_pden(s,t,p,0);
sgth_fit=sw_pden(polyval(P,t),t,p,0);

subplot(224)
hraw = plot(sgth,p);
hold on
hfit = plot(sgth_fit,p,'linewidth',2)
axis ij
grid on
xlabel('sgth')
ylabel('P')
legend([hraw hfit],'raw','fit','location','best')


%% For patches where R^2>0.5, compute s from T using fit, then compute N2, dT/dz etc

igc = find(patches.R2>0.5);

patches_good      = struct();
patches_good.p1   = patches.p1(igc) ;
patches_good.p2   = patches.p2(igc) ;
patches_good.chi  = patches.chi(igc) ;
patches_good.eps  = patches.eps(igc) ;
patches_good.cnum = patches.cnum(igc) ;
patches_good.R2   = patches.R2(igc) ;
patches_good.n2_line   = patches.n2_line(igc) ;
patches_good.dtdz_line   = patches.dtdz_line(igc) ;

patches_good.n2_fit = nan*ones(size(patches_good.p1));
cnum_loaded=-5;
hb = waitbar(0)
for ip=1:length(igc)
    waitbar(ip/length(igc),hb)
    if cnum_loaded~=patches_good.cnum(ip)
        % load raw chameleon cast
        clear cal cal2 head
        cal = load_cal_tiwe(cnum) ;
        cnum_loaded = cnum;
    else
    end
    
    % get data in patch
    clear idz t p sraw P S s ptmp sgth
    idz = isin(cal.P,[patches_good.p1(ip) patches_good.p2(ip) ]);
    
    if length(idz)>10
        
        t = cal.T1(idz)  ;
        p = cal.P(idz)   ;
        sraw = cal.SAL(idz) ;
        
        [P,S] = polyfit(t,sraw,1);
        s = polyval(P,t);
        %
        ptmp=sw_ptmp(s,t,p,0);
        sgth=sw_pden(s,t,p,0);
        
        % sorth pot. dens.
        clear sgth_sort I
        [sgth_sort , I]=sort(sgth,1,'ascend');
        
        %~~ compute drho/dz and N^2
        % fit a line to sgth
        clear P1
        P1 = polyfit(p,sgth,1);
        
        % calculate N^2 from this fit
        clear drhodz
        drhodz = -P1(1);
        patches_good.n2_fit(ip) = -9.81/nanmean(sgth)*drhodz;
        
    end
    
end
delete(hb)

%% compare N^2 computed from fit to that from raw data

figure(1) ; clf
%loglog(patches_good.n2_line,patches_good.n2_fit,'.')
%scatter( log10(patches_good.n2_line), log10(patches_good.n2_fit),'filled','MarkerFaceAlpha',0.1)
histogram2( log10(patches_good.n2_line), real(log10(patches_good.n2_fit)),50,'DisplayStyle','tile')
grid on
xlim([-5.5 -3])
ylim([-5.5 -3])
hold on
xvec=linspace(-5.5,-3,100);
plot(xvec,xvec,'k--')
xlabel('N^2 line')
ylabel('N^2 fit')



%%
figure(2);clf
histogram(real(log10(patches_good.n2_fit./patches_good.n2_line)))
%%

figure(2);clf
histogram(log10(patches_good.n2_line),'Normalization','pdf')
hold on
histogram(real(log10(patches_good.n2_fit)),'Normalization','pdf')

%%

gam1 = ComputeGamma(patches_good.n2_line, patches_good.dtdz_line, patches_good.chi, patches_good.eps);
gam2 = ComputeGamma(patches_good.n2_fit , patches_good.dtdz_line, patches_good.chi, patches_good.eps);

figure(1);clf
h1=histogram(log10(gam1),'Normalization','pdf','EdgeColor','none')
hold on
h2=histogram(real(log10(gam2)),'Normalization','pdf','EdgeColor','none')
grid on
xlim([-3 2])

freqline(log10(0.2))
xlabel('log_{10}\gamma')

% plot Bill's TIWE results too
load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/tiwe_patch_gamma/data/events_TIWE.mat')
gamA=ComputeGamma(A.N2,A.tgrad,A.chi,A.eps);
hold on
h3=histogram(real(log10(gamA)),'Normalization','pdf','EdgeColor','none')

legend([h1 h2 h3 ],'raw','fit','Bill','location','best')

%%


