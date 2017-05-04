function [chipod, cham]=discard_convection_eq14(chipod,cham)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14_patch_gamma/data/EQ14_mldepths.mat')

for i=1:length(cham.cnum)
    try
        id = find(cham.cnum(i)==zml_cnum) ;
        izbad = find(cham.P<zml(id));
        cham.eps(izbad,i)=nan;
    end
end

for i=1:length(chipod.cnum)
    try
        id = find(chipod.cnum(i)==zml_cnum) ;
        izbad = find(chipod.P<zml(id));
        chipod.eps(izbad,i)=nan;
    end
end

