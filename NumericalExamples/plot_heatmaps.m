%% Plot the HeatMaps
% This scripts plots the HeatMap generated with the heatmaptest script.
% The script is made to run in batch mode, hence the analysis of the
% results is done here.

clear; clc; close all;

% Load the data from file
load('heatmapresults.mat')

% Load the list of test matrices
QP_problems_path = "../testmatrices/";
d = dir(fullfile(QP_problems_path,'*.mat')); 

% Data:
perc = [0.1,0.2,0.3,0.4,0.5];
betav = [0.8,0.6,0.4,0.2,1];


nnzA = zeros(length(d),1);
JA = zeros(length(perc),length(betav));
for i = 1:length(d)
    load(fullfile(QP_problems_path,d(i).name));
    A = Problem.A;
    nnzA(i) = nnz(A);
    for j = 1:length(betav)
        JA(i,j) = betav(j)*norm(A,"fro")^2 + (1-betav(j))*norm(A,1);
    end
end


% Performances in term of the objective function
figure('Position',[958 838 1220 287],'Name',"Objective Function")
for i=1:5
    subplot(1,5,i)
    heatmap(log10(Pert(:,:,i)./JA), ...
        "Colormap",turbo(30),...
        "XData",perc, ...
        "XLabel","Fraction of top nodes averaged", ...
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
        "XLabel","Fraction of top nodes averaged", ...
        "Title",sprintf("\\beta = %1.1f",betav(i)));
end

% Number of nonzero entries in the perturbation
figure('Position',[958 838 1220 287],'Name',"nnz(\Delta)")
for i=1:5
    subplot(1,5,i)
    heatmap(reshape(NNZ(:,i,:)./nnzA,15,5,1), ...
        "Colormap",turbo(40),...
        "XData",betav, ...
        "XLabel","\beta", ...
        "Title",[sprintf("x = %1.1f %%",perc(i)*100)]);
end