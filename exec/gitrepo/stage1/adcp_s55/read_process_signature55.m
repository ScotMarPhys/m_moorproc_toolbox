%%  Code for the quality control of Signature 55 ADCP

% K Burmeister, S Jones 12/2025


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PARAMETER PRAEMBLE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% Need to be updated for:
% - directories and file names
% - add additional deployment periods (add moorX = 'filename')
% - turn off/on checkplots
% - edit start/end of total time period
% - data version
% - set depth of shallowest instrument (idepth)
% - preamble for despiking
% - if you add a deployment period, you need to edit steps 1-3
%---------
close all
clear

addpath(genpath('../functions'));

% TO UPDATE <------------------------
% in- and output directories
pc_name = getenv('COMPUTERNAME');
if strcmp(pc_name,'SA07KB-3JN9YY2');
    basedir = 'C:\Users\sa07kb\OneDrive - SAMS\data\data_OSNAP\German_53N_Array\';
    pathgit = 'C:\Users\sa07kb\Projects\Oi\';
elseif strcmp(pc_name,'SA01SJ-G9WC2J3')
    dataindir = 'E:\OSNAP\RHADCP\DY181\S200044A012_RHAD2_JC238\conversion2\';
    pathgit = 'D:\Work_computer_sync\OSNAP_postdoc\Python\m_moorproc_toolbox\';    
    figureoutdir = ['D:\Work_computer_sync\OSNAP_postdoc\Mooring\RHADCP\plots\'];
    addpath(genpath('D:\Work_computer_sync\MATLAB_functions')); % General functions
else
    error('Please add your path above')
end

filename = 'S200044A012_RHAD2_JC238';
%% Load
load([dataindir filename '.mat']);


%% Build post-processing flag array

% 0: QC_NOT_EVALUATED
% 1: QC_GOOD
% 2: QC_UNKNOWN
% 3: QC_PROBABLY_BAD
% 4: QC_BAD
% 5: QC_CHANGED
% 6: QC_UNSAMPLED
% 7: QC_INTERPOLATED
% 8: ?
% 9: QC_MISSING

QC_vel = 0*double(Data.Average_VelEast);

% NOTE: Do we want the ability to interactively flag and zoom on these early plots?
% E.g. It's a good opportunity to remove deployment and recovery based on
% pressure etc.  A bit of automatic flagging for demo in Fig. 1

%% STAGE 1.  Nortek suggested quality control steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timemin = Data.Average_Time(1);
timemax = Data.Average_Time(end);

%% 1. Tilt / pitch

figure(1); 

title([filename ' pitch'],'fontsize',14,'Interpreter','none');
hold on; grid on;
% plot data
plot(Data.Average_Time,Data.Average_Pitch,'.');

datetick('keepticks','keeplimits');

% Suggested quality thresholds
line([timemin timemax],[0 0],'color','k');
line([timemin timemax],[10 10],'color','g');
line([timemin timemax],[-10 -10],'color','g');
line([timemin timemax],[30 30],'color','r');
line([timemin timemax],[-30 -30],'color','r');

text(timemin+50,12,'10 < 30: Post processing possible','color','g');
text(timemin+50,-12,'-10 > -30: Post processing possible','color','g');
text(timemin+50,32,'>30: Profiles probably bad','color','r');
text(timemin+50,-32,'<-30: Profiles probably bad','color','r');

xlim([timemin timemax]);
ylim([-40 40]);

xlabel('Date');
ylabel('Pitch (^o)');

% FLAG EXAMPLE: Automatically flag any bad timesteps
badind = find(Data.Average_Pitch > 30 | Data.Average_Pitch < -30);
plot(Data.Average_Time(badind),Data.Average_Pitch(badind),'or');
QC_vel(badind,:) = 4;
disp([num2str(length(badind)) ' timesteps flagged bad (4) due to excessive pitch']);

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f1_pitch_QC']);

%% heading %%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
title([filename ' heading'],'fontsize',14,'Interpreter','none');
hold on; grid on;

% Plot data
plot(Data.Average_Time,Data.Average_Heading,'g.');
plot(Data.Average_Time,Data.Average_Heading,'g');

datetick('keepticks','keeplimits');
xlim([timemin timemax]);
ylim([-0 360]);
xlabel('Date');
ylabel('Heading (^o)');

% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f2_heading_QC']);


%% Pressure %%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
title([filename ' pressure'],'fontsize',14,'Interpreter','none');
hold on; grid on;

% Plot data
plot(Data.Average_Time,Data.Average_Pressure,'k.');

datetick('keepticks','keeplimits');
xlim([timemin timemax]);
% ylim([-0 360]);
xlabel('Date');
ylabel('Pressure (db)');
set(gca,'YDir','reverse');


% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12]*1.5)
print('-dpng',[figureoutdir filename '_f3_pressure_QC']);




