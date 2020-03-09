function [D] = distanceMat(Z, A)
%distanceMat - calculate Euclidean distance between each pair embeddings
% Syntax: [D] = distanceMat(Z, A)
%
% Inputs:
%   Z - nodeNum*m, z_i is the embeding vector of v_i
%   A - nodeNum*nodeNum, the adjacency matrix of the network
% 
% Outputs:
%   D - nodeNum*nodeNum, d_ij is the Euclidean distance of z_i and z_j
%
% Author:  Y. Dong
% Created: Jun 21, 2019

% Variables definition
nodeNum = size(A, 1);
D = zeros(nodeNum, nodeNum);

% calculate the distance
for i = 1:nodeNum
    for j = 1:nodeNum
        if A(i,j) == 1
            D(i,j) = sqrt(sum((Z(i,:)-Z(j,:)).^2));
        end
    end
end
%D(D==0) = 1000*max(max(D));
D_self = max(D);
D = D - diag(diag(D)) + diag(D_self);
end