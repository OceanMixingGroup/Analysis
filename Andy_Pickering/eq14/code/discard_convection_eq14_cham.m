function avg = discard_convection_eq14_cham(avg,cnum)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Discard chi and epsilon data from mixed layer, identified in
% Identify_ML_eq14.m
%
%
% INPUT
% - avg  : binned Chameleon profile
% - cnum : castnumber
%
% OUTPUT
% - avg : Chameleon profile w/ ML data NaNd out
%
%~~~~~~~~~~~~
% A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

load('/Users/Andy/Cruises_Research/Analysis/Andy_Pickering/eq14/data/EQ14_mldepths.mat')

try
    id = find( cnum==zml_cnum ) ;
    izbad = find(avg.P<zml(id));
    avg.CHI(izbad)=nan;
    avg.EPSILON(izbad)=nan;
end

%disp([num2str(length(izbad)) ' points removed']);
