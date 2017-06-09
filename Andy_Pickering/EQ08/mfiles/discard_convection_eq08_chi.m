function avg = discard_convection_eq08_chi(avg,cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq08_patch_gamma/data/eq08_mldepths.mat')

try
    id = find( cnum==zml_cnum ) ;
    izbad = find(avg.P<zml(id));
    avg.chi1(izbad)=nan;
    avg.eps1(izbad)=nan;
end

%disp([num2str(length(izbad)) ' points removed']);

