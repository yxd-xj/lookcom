%-----------------------------------------------------------
% Title:   Synthetic Expriment Plot
% Author:  Y. Dong
% Created: Jul 21, 2019

% Description: Draw adjacency matrix and the optimal network structure
%-----------------------------------------------------------

% Variables Definition
nodeNum = size(A,1);
x_lable = cell(nodeNum,1);
x_lable(:,1) = {''};
y_lable = cell(nodeNum,1);
y_lable(:,1) = {''};

% plot the adjacency matrix
subplot(1,2,1);
h1 = heatmap(A, 'XLabel','','YLabel','');
h1.XDisplayLabels = x_lable;
h1.YDisplayLabels = y_lable;
h1.GridVisible = 'off';
title("Original Network Structure A");

% plot the learnged optimal network structure
subplot(1,2,2);
h2 = heatmap(W,'XLabel','','YLabel','');
colorbar;
h2.XDisplayLabels = x_lable;
h2.YDisplayLabels = y_lable;
h2.GridVisible = 'off';
title("Optimal Network Structure W");