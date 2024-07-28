%% Plot the HeatMaps
% This scripts plots the HeatMap generated with the heatmaptest script.
% The script is made to run in batch mode, hence the analysis of the
% results is done here.

clear; clc; close all;

% Load the data from file
load('heatmapresults.mat')

perc = [0.1,0.2,0.3,0.4,0.5];
betav = [1,0.8,0.6,0.4,0.2];

% Performances in term of the objective function
figure('Position',[958 838 1220 287],'Name',"Objective Function")
for i=1:5
    subplot(1,5,i)
    heatmap(log10(Pert(:,:,i)), ...
        "Colormap",turbo(30),...
        "XData",perc, ...
        "XLabel","Percentage of top nodes averaged", ...
        "Title",sprintf("\\beta = %1.1f",betav(i)));
end

% Kendall's tau rank correlation score between the desired and
% obtained measure
figure('Position',[958 838 1220 287],'Name',"Kendall's \tau")
for i=1:5
    subplot(1,5,i)
    heatmap(KTAU(:,:,i), ...
        "Colormap",turbo(10),...
        "XData",perc, ...
        "XLabel","Percentage of top nodes averaged", ...
        "Title",sprintf("\\beta = %1.1f",betav(i)));
end

% Number of nonzero entries in the perturbation
figure('Position',[958 838 1220 287],'Name',"nnz(\Delta)")
for i=1:5
    subplot(1,5,i)
    heatmap(reshape(log10(NNZ(:,i,:)),15,5,1), ...
        "Colormap",turbo(40),...
        "XData",betav, ...
        "XLabel","\beta", ...
        "Title",[sprintf("x = %1.1f %%",perc(i)*100)]);
end