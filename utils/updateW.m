function [W, neighborsNum] = updateW(Z, F, A, LAMBDA, ALPHA)
%updateW - update W by solving the following problem:
%   min  sum_{i=1}^n ||z_i - \sum_{j=1}^N w_ij z_j||_2 + sum_{i=1}^n alpha_i ||w_i - p_i||_2^2
%        + LAMBDA/2 sum_{i,j=1}^n w_ij*||f_i - f_j||_2^2
%   s.t. 1'w_i = 1, w_i >= 0
%
% Syntax: [W, neighborsNum] = updateW(Z, F, LAMBDA, ALPHA)
%
% Inputs:
%   Z - nodeNum*m, representation matrix of n nodes
%   F - nodeNum*clusterNum, matrix for low rank constraint
%   A - nodeNum*nodeNum, the adjacency matrix of the network
%   LAMBDA - float, parameter for low rank constraint
%   ALPHA - float, parameter to control the similarity bettween the learned structure and the original structure
%  
% Outputs:
%   W - nodeNum*nodeNum, optimal community structure of the network
%   neighborsNum - nodeNum*1, number of the selected neighbors of nodes
%
% Author:  Y. Dong
% Created: Jun 24, 2019

% Variables definition
nodeNum = size(A, 1);
W = zeros(nodeNum, nodeNum);
neighborsNum = zeros(nodeNum, 1);

% construct the distance matrix
D_z = distanceMat(Z, A);
D_f = distanceMat(F, A);
D = D_z.^2 + 0.5*LAMBDA*(D_f.^2);

% update w_i for each node
for i = 1:nodeNum
    % find the neighbors of node v_i
    neighbor_idx = find(A(i,:)>0);
    d_neighbor = D(i, neighbor_idx);
    % sort neighbors by D
    [d_neighbor_sorted, sorted_idx] = sort(d_neighbor);
    d_neighbor_sorted = [d_neighbor_sorted 10*max(d_neighbor_sorted)];
    % get p_i according to d_nieghbor_sorted
    p = A(i,:)./sum(A(i,:));
    p_neighbor_sorted = p(neighbor_idx(sorted_idx));
    p_neighbor_sorted = [p_neighbor_sorted 0];
    % define variables
    theta = d_neighbor_sorted(1) - 2*ALPHA*p_neighbor_sorted(1) + 2*ALPHA;
    neighbors = 0;
    neighbors_num = length(d_neighbor);
    sum_d = 0;
    sum_p = 0;
    % calculate theta and k^*
    while ((2*ALPHA*p_neighbor_sorted(neighbors+1)+theta) > d_neighbor_sorted(neighbors+1)) && (neighbors < neighbors_num)
        neighbors = neighbors + 1;
        sum_d = sum_d + d_neighbor_sorted(neighbors);
        sum_p = sum_p + p_neighbor_sorted(neighbors);
        theta = 2*ALPHA/neighbors - 2*ALPHA/neighbors*sum_p + 1/neighbors*sum_d;
    end
    % update w_i
    neighborsNum(i) = neighbors;
    W(i, neighbor_idx(sorted_idx(1:neighbors))) = p(neighbor_idx(sorted_idx(1:neighbors))) + 1/(2*ALPHA)*theta - 1/(2*ALPHA)*D(i, neighbor_idx(sorted_idx(1:neighbors)));
    W(i,:) = W(i,:)./sum(W(i,:));
end

end