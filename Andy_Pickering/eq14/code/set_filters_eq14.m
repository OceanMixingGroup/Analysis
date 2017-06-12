% set_filters.m
% A suppliment file for average_data_gen.
% Here you can set filters that will be used to calculate
% EPSILON in calc_epsilon_filt_gen.m routine.
% In addition high and low frequency cutoff and high frequency 
% cutoff coefficient are defined here (see documentation to 
% calc_epsilon_filt_gen.m).
% Choice of filters and coefficients depends on cast number
% (it could be different for different tows)
% Note: variable CAST (cast number) should be defined as GLOBAL
% in main script and in average_data_gen.m routine

%%%%%% notes by Sally (March 2015)
%%% Filters are set to get rid of any narrowband strumming that may be
%%% occurring as the chameleon is falling such as the EUC flowing past the
%%% line or the crash ring shedding eddies. It's important that these are
%%% filtered out so they don't get into the final Epsilon signal. How to
%%% know what's right? Usually want the smallest values of epsilon in the
%%% averaged profile. Looking at individual spectra at different depths for
%%% different instrument configurations is really the only way to tell what
%%% is BEST.



%$$$$$$$$ DO NOT CHANGE THIS SECTION!! $$$$$$$$$$$$$$$$$$$$$
global q
cast=q.script.num;
clear filters;
% SET DEFAULTS
filters={};
k_start=2; % [cpm]
k_stop=90; % [cpm]
kstop_coef=0.5;
% eps=eval(['avg.EPSILON' prb '(n)']);
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%##### IN THIS SECTION YOU SET FILTERS AND COEFFICIENTS ####
%##### AS FUNCTION OF CAST NUMBER ##########################
% Filters set according to mean of first 4 points of shear spectra
% Filters could be different for diferent probes (variable prb)
% and could depend on ADP (variable cal.FLAG, cal.FLAG(:)=1 when 
% ADP is ON)
k=mean(log10(ss(1:4)));
% if cast>=78
%     if k<-5
%         filters={'l7'};
%     elseif k<-4.8
%         filters={'l10'};
%     elseif k<-4.6
%         filters={'l12'};
%     elseif k<-4.4
%         filters={'l15'};
%     elseif k<-4.0
%         filters={'l20'};
%     elseif k<-3.7
%         filters={'l22'};
%     else
%         filters={'l40'};
%     end
% else

%     if k<-5
%         filters={'l6'};
%     elseif k<-4.8
%         filters={'l8'};
%     elseif k<-4.6
%         filters={'l10'};
%     elseif k<-4.4
%         filters={'l12'};
%     elseif k<-4.2
%         filters={'l15'};
%     elseif k<-4.0
%         filters={'l18'};
%     elseif k<-3.7
%         filters={'l20'};
%     else
%         filters={'l40'};
%     end
    
    
%     if k<-5
%         filters={'l7'};
%     elseif k<-4.8
%         filters={'l8'};
%     elseif k<-4.6
%         filters={'l10'};
%     elseif k<-4.4
%         filters={'l12'};
%     elseif k<-4.0
%         filters={'l20'};
%     elseif k<-3.7
%         filters={'l22'};
%     else
%         filters={'l40'};
%     end
%     
%     
    
    

%%% these are the filters from EQ08. Due to the undercurrent, having
%%% depth-dependant filters makes sense. (Dynamo filters were not depth
%%% dependant.) Sasha would change the filters throughout the cruise based
%%% on the fall speed and configuration of the instrument. However, for
%%% now, I will assume the setup was similar to that for EQ08 and use these
%%% filters for the entire EQ14 cruise.
if avg.P(n)<80
	if k<-4.4
        filters={'l13'};
    elseif k<-4.1
        filters={'l18'};
%             filters={'l15'};
%              filters={'l10'};
   elseif k<-3.5
        filters={'l22'};
%             filters={'l15'};
    else
        filters={'l28'};
%             filters={'l15'};
%             filters={'l10'};
    end
elseif avg.P(n)>=80 && avg.P(n)<160
    if k<-4.4
        filters={'l7'};
    elseif k<-4.1
        filters={'l10'};
    elseif k<-3.5
        filters={'l13'};
    else
        filters={'l20'};
%             filters={'l15'};
%             filters={'l10'};
    end
else
    if k<-4.4
        filters={'l4'};
    elseif k<-4.1
        filters={'l8'};
    elseif k<-3.5
        filters={'l9'};
    else
        filters={'l20'};
%             filters={'l15'};
%             filters={'l10'};
    end
end
    