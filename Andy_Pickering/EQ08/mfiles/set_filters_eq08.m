% set_filters.m
% THESE ARE THE FILTERS FOR CR04!!! THEY ARE NOT CORRECTED FOR CR05!!!
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
% Filters set according to EPSILON estimation (variable eps)
% EPSILON estimated in average_data_gen script before applying filters
% Filters could be different for diferent probes (variable prb)
% and could depend on ADP (variable cal.FLAG, cal.FLAG(:)=1 when 
% ADP is ON)
k=mean(log10(ss(1:4)));
head_index_num=eval(['head.sensor_index.S' prb]);
if q.script.num<1000
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
else
    if avg.P(n)<80
        if k<-4.4
            filters={'l13'};
        elseif k<-4.1
            filters={'l18'};
%             filters={'l15'};
%             filters={'l10'};
        elseif k<-3.5
            filters={'l24'};
%             filters={'l15'};
%             filters={'l10'};
        else
            filters={'l28'};
%             filters={'l15'};
%             filters={'l10'};
        end
    else
        if k<-4.4
            filters={'l4'};
        elseif k<-4.1
            filters={'l5'};
        elseif k<-3.5
            filters={'l10'};
        else
            filters={'l18'};
%             filters={'l15'};
%             filters={'l10'};
        end
    end
end
%#########################################################    
