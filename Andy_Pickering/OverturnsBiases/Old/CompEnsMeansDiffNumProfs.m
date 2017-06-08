%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompEnsMeansDiffNumProfs.m
%
% Goal: Previously I have been computing ensemble mean from resampled data
% using 100 profiles shifted by 2 min each, for all diferent speeds. But I
% think what I really should be doing is using the # of shifted profiles
% that goes through the time taken for 1 profile (ie if the profile takes 1
% hour, should use 30 profiles shifted by 2 mins.
%
% In this code, I will compute hte ensemble means both ways and see if it
% makes a big difference.
%
% Doesn't seem to make a big difference. Turns out 100 is close to a
% multiple of T_prof / w_samp.
%
% 16 Jan. 2015
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

whmoor=3
testnum=4

load(fullfile(cd,'Data',['Tchain' num2str(whmoor) '_Test' num2str(testnum) '_AllCases']))

zrange=range(REsamp.z)
T=zrange/REsamp.w_samp/60 % mins
np=round(T/2)

figure(1);clf

em=squeeze(nanmean(REsamp.eps,2));

semilogx(nanmean(em,2),REsamp.z)
axis ij
hold on
semilogx(nanmean(em(:,1:np),2),REsamp.z)
xlim([1e-8 1e-5])
%%