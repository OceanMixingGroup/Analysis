function [out]=compute_Tz_N2_for_patch(p1,p2,p,t,s,ptmp,sgth,alpha,Lt)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% function to compute dT/dz and N^2 etc for a single patch
%
% INPUTS
% p1,p2,p,t,s,ptmp,sgth, alpha, Lt
%
% OUTPUTS
%
% dtdz_range
% dtdz_line
% dtdz_bulk
%
% n2_range
% n2_line
% n2_bulk
% n4
%
%-------------------
% 2/28/17 A. Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% get indices of data in patch
clear iz
iz = isin(p,[ p1 p2 ]) ;

if length(iz)>10
    
    clear t_ot s_ot p_ot ptmp_ot sgth_ot
    t_ot = t(iz);
    s_ot = s(iz);
    p_ot = p(iz);
    ptmp_ot = ptmp(iz);
    sgth_ot = sgth(iz);
    
    % sort temp
    clear t_sort I
    [t_sort , I] = sort(t_ot,1,'descend');
    
    % sort potential temp
    clear t_pot t_pot_sort
    [ptmp_sort , Iptmp]=sort(ptmp_ot,1,'descend');
           
    % fit a line to pot. temp to get slope
    clear P dtdz_line
    P = polyfit(p_ot,ptmp_sort,1);
    dtdz_line=-P(1);
    
    clear P
    
    %~~ 'bulk gradient' method from Smyth et al 2001
    % essentially = rms T (btw sorted/unsorted) /  thorpe scale ?
    %    t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
    clear t_rms dtdz_bulk
    t_rms = sqrt( nanmean(( ptmp_ot - ptmp_sort ).^2) );
    dtdz_bulk = t_rms / Lt ;
    
    %~~ Now do similar for density / N^2
    
    % sorth pot. dens.
    [sgth_sort , I]=sort(sgth_ot,1,'ascend');
        
    % fit a line to sgth
    clear P1
    P1 = polyfit(p_ot,sgth_sort,1);
    
    % calculate N^2 from this fit
    clear drhodz n2_2 drho dz n2_3
    drhodz = -P1(1);
    drhodz_line = drhodz;
    n2_line = -9.81/nanmean(sgth)*drhodz;
    
    % compute bulk N^2 (need to fit 
    n2_bulk= -9.81 / nanmean(sgth) * alpha * t_rms / Lt    ;    
    
    % compute N^2 w/ sw_bfrq
    clear n2
    n2 = sw_bfrq(s_ot(I),t_ot(I),p_ot,0.03);
    n4=nanmean(n2);    

    % Organize results in 'out' structure
    out.dtdz_line  = dtdz_line;
    out.dtdz_bulk  = dtdz_bulk;
    
    out.n2_line  = n2_line;
    out.n2_bulk  = n2_bulk;
    out.n4 = n4;
    
end

return

%%