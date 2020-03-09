function [toTerminate, obj] = terminate(Z, F, W, A, ALPHA, LAMBDA, obj_pre)
%terminate - Determine whether to terminate the update
% Syntax: [toTerminate, obj, LAMBDA] = terminate(Z, F, W, A, ALPHA, LAMBDA_init, LAMBDA_cur, obj_pre, eigvals_F, clusterNum)
%
% Inputs:
%   Z - nodeNum*m, representation matrix of n nodes
%   F - nodeNum*clusterNum, matrix for low rank constraint
%   W - nodeNum*nodeNum, optimal community structure of the network
%   A - nodeNum*nodeNum, the adjacency matrix of the network
%   ALPHA - float, parameter to control the similarity bettween the learned structure and the original structure
%   LAMBDA - float, parameter for low rank constraint
%   obj_pre - float, the last objective function value
%  
% Outputs:
%   toTerminate - bool, flag to indicate whether to terminate the update
%   obj - float, the objective function value
% 
% Author:  Y. Dong
% Created: 12 Feb, 2020

% Variables definition
TOL = 5e-3;
toTerminate = false;
nodeNum = size(A, 1);
P = zeros(nodeNum, nodeNum);
for i = 1:nodeNum
    degree_i = sum(A(i,:));
    P(i,:) = A(i,:)./degree_i;
end

% Calculate the objective function value
obj = sum(sum(W.*(distanceMat(Z,A).^2))) + ALPHA*(norm(W-P, 'fro'))^2 + LAMBDA*(sum(sum(W.*(distanceMat(F,A).^2))));

% determine whether to terminate the update
err = abs(obj - obj_pre);
if (err < TOL) || (obj_pre < obj)
    toTerminate = true;
end

end