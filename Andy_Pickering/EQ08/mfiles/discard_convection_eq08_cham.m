function avg = discard_convection_eq08_cham(avg,cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq08/data/eq08_mldepths.mat')

try
    id = find( cnum==zml_cnum ) ;
    izbad = find(avg.P<zml(id));
    avg.CHI(izbad)=nan;
    avg.EPSILON(izbad)=nan;
end

%disp([num2str(length(izbad)) ' points removed']);

