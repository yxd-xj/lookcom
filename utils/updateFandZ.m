function [F, eigvals] = updateFandZ(W, c)
%updateF - Update F by solving the following problem:
%       min_F   Tr(F'L_wF)
%       s.t.    F \in R^{n*c}, F'F = I_c
%
% Syntax: [F, eigvals] = updateF(W, c)
%
% Inputs:
%   W - nodeNum*nodeNum, the optimal network structure for cd
%   c - int, number of communities on the network
%
% Outputs:
%   F - nodeNum*c, matrix for low rank constraint
%   eigvals - c*1, c smallest eigenvalues of L_W
%
% Author:  Y.X. Dong
% E-mail:  dyx_xjtu@163.com
% Created: Jun 21, 2019

% Variables definition
ZR = 1e-10;
nodeNum = size(W, 1);

% console
isWnotSym = true;

% Construct the Laplacian matrix
if isWnotSym
    A = max(W', W);
    % A = (W' + W)./2;
else
    A = W;
end
D = diag(sum(A));
L = D - A;
L_rw = diag(1./(diag(D)))*L;

% Get c smallest eigenvalues and eigenvectors
[F, eigvals] = eigs(L_rw, c+1, -1);
eigvals = diag(eigvals);
[eigvals, idx] = sort(eigvals);
eigvals(abs(eigvals)<ZR) = eps;
eigvals = eigvals(1:c+1);
F = F(:, idx(1:c));
% F = F*diag(sqrt(1./diag(F'*F)));
% F = normalize(F, 2, 'norm');

end