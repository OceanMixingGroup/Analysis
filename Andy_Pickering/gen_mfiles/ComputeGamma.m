function gam=ComputeGamma(n2,dtdz,chi,eps)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute gamma (mixing efficiency)
% 
%~~~~~~~~~~~~~~
% 10/27/16 - APickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
gam =  n2 .* chi ./2 ./ eps ./ (dtdz.^2);

return
%%
