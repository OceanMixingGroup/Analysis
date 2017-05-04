function avg = discard_convection_eq14_chi(avg,cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat')

try
    id = find( cnum==zml_cnum ) ;
    izbad = find(avg.P<zml(id));
    avg.eps1(izbad)=nan;
end

%disp([num2str(length(izbad)) ' points removed']);

