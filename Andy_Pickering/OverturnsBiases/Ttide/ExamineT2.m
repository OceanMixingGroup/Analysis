%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ExamineT2.m
%
% hi Andy - I made a datafile for you from ttide.  assume S is constant; 
% T is in-situ, distances are heights above the bottom and the bottom depth
% is something around 1900 m?  depths are not correct, but fine for your 
% tests right now I suspect
%
% 28 Feb 2015
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Dropbox/ttide_jn/T2_short_timeseries.mat')

%

id=find(dat.T>3.5);
dat.T(id)=nan;

tm=nanmean(dat.T,2);
zc=dat.depth(1):-10:dat.depth(end);
tmi=interp1(dat.depth,tm,zc);

figure(1);clf
ezpc(datenum2yday(dat.time),dat.depth,dat.T)
hold on
contour(datenum2yday(dat.time),dat.depth,dat.T,tmi(1:7:end),'k')
colorbar
axis xy
xlabel('Yearday')
ylabel('Height Above Bottom')

%%
z=1900-dat.depth;
figure(2);clf
ezpc(datenum2yday(dat.time),z,dat.T)
hold on
contour(datenum2yday(dat.time),z,dat.T,tmi(1:7:end),'k')
colorbar
axis ij
xlabel('Yearday')
ylabel('Depth')
caxis([2.1 3.3])

%%

fname=fullfile('/Users/Andy/Cruises_Research/IWISE/Analysis/S9/Dissipation/OverturnsBiases','NotesOverturnBiases',['T2_Temp'])
addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/exportfig/
export_fig(fname,'-png')

%fname=
%save2pdf(fname)

%%