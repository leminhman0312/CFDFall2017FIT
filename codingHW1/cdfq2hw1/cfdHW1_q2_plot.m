%% Import data from text file.
% Script for importing data from the following text file:
%
%    /home/maxle/Dropbox/FIT GRAD MASTERS/FALL 2017/CFD_FALL_17/codingHW1/cdfq2hw1/CFD_hw1q2_ResultsPlottingMatlab.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/10/01 21:24:58

%% Initialize variables.
filename = '/home/maxle/Dropbox/FIT GRAD MASTERS/FALL 2017/CFD_FALL_17/codingHW1/cdfq2hw1/CFD_hw1q2_ResultsPlottingMatlab.txt';
delimiter = '\t';

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);



%% Close the text file.
fclose(fileID);



%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
CFDhw1q2ResultsPlottingMatlab = [dataArray{1:end-1}];

%%Prepare X and Y variables to plot

x = 0:0.1:3.5;
y = 0:0.1:3.5;
T = CFDhw1q2ResultsPlottingMatlab;

contourf(x,y,T,20,'LineColor','none');
xlabel('X [ft]','fontsize',18);
ylabel('Y [ft]','fontsize',18);
box on 
title('Solving $$\frac{\partial T}{\partial t} = {\alpha}\bigg[\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}\bigg]$$ at t = 0.4 hours on a 3.5ft x 3.5ft square,chrome steel bar','interpreter','latex','fontsize',18);
colormap(jet)
shading flat
tickvector = [0,50,100,150,200];
colorbar('Direction','reverse','Ticks',tickvector);
set(gca,'fontsize',20)
k = waitforbuttonpress




%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
