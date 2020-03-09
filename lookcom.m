function [Z, W, F, neighborsNum, objs] = lookcom(A, clusterNum, ALPHA, LAMBDA, m)
%lookcom - learning network embedding for community detection
%
% Syntax: [Z, W, F, neighborsNum, objs] = cdOptNetEm(A, clusterNum, ALPHA, LAMBDA, m)
%
% Inputs:
%   A - nodeNum*nodeNum, the adjacency matrix of the network
%   clusterNum - int, number of communities on the network
%   ALPHA - float, parameter to control the similarity bettween the learned structure and the original structure
%   LAMBDA - float, parameter for low rank constraint
%   m - int, dimension of node representation
%
% Outputs:
%   Z - nodeNum*m, representation matrix of n nodes
%   W - nodeNum*nodeNum, optimal community structure of the network
%   F - nodeNum*clusterNum, matrix for low rank constraint
%   neighborsNum - nodeNum*1, number of the selected neighbors of nodes
%   objs - MAXITERS*1, the objective function values in optimazation
%
% Author:  Y. Dong
% Created: Jun 20, 2019

% Settings
ZR = 10e-11;
TOL = 5e-3;
MAXITERS = 50;
LAMBDA_init = LAMBDA;
nodeNum = size(A, 1);
objs = zeros(MAXITERS, 1);

% Console
isInitW = true;
isWnotClear = true;

% initialize
if isInitW
    P = zeros(nodeNum, nodeNum);
    for i = 1:nodeNum
        degree_i = sum(A(i,:));
        P(i,:) = A(i,:)./degree_i;
    end
    W = P;
else
    W = A;
end
W_org = W;
Z = zeros(nodeNum, m);
F = zeros(nodeNum, clusterNum);

% the main algorithm
for iter = 1:MAXITERS
    fprintf("Iteration %d start...\n", iter);
    % update F
    F_pre = F;
    [F, eigvals_F] = updateFandZ(W, clusterNum);
    if iter > 1
        [isWnotClear, F_eigvals_sum_c_plus] = updateStatus(eigvals_F, clusterNum, F_eigvals_sum_c_plus);
        if isWnotClear
            [toTerminate, obj] = terminate(Z, F, W_org, A, ALPHA, LAMBDA_init, obj);
            if toTerminate
                F = F_pre;
            end
        end
    else
        obj = 0;
        F_eigvals_sum_c = sum(eigvals_F(1:clusterNum));
        F_eigvals_sum_c_plus = sum(eigvals_F(1:clusterNum+1));
        fprintf("f_c = %f; f_c_plus = %f \n", F_eigvals_sum_c, F_eigvals_sum_c_plus)
    end
    % update Z
    Z_old = Z;
    [Z, eigvals_Z] = updateFandZ(W, m);
    if isWnotClear
        [toTerminate, obj] = terminate(Z, F, W_org, A, ALPHA, LAMBDA_init, obj);
    else
        F = F_pre;
        [toTerminate, obj] = terminate(Z, F, W_org, A, ALPHA, LAMBDA_init, obj);
        W = max(W, W');
        break;
    end
    % update W
    W_old = W;
    [W, neighborsNum] = updateW(Z, F, A, LAMBDA, ALPHA);
    W_org = W;
    [toTerminate, obj] = terminate(Z, F, W_org, A, ALPHA, LAMBDA_init, obj);
    objs(iter) = obj;
    if toTerminate
        W = max(W, W');
        break;
    end
end
end