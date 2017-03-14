function R2 = compute_R2(x,y,P)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute R^2 for a linear regression fit of y~x
%
% P is result of polyfit(x,y)
%
%
% 3/14/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
yfit = polyval(P,x);
yresid = y-yfit;
SSresid=sum(yresid.^2);
SStotal=(length(y)-1) *var(y);
R2 = 1 - SSresid/SStotal;
%%