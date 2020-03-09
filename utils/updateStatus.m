function [isWnotClear, F_eigvals_sum_c_plus] = updateStatus(eigvals_F, clusterNum, F_eigvals_sum_c_plus_pre)
%myFun - deterimate whether to update W and update LAMBDA
% Syntax: [update_w, LAMBDA, F_eigvals_c] = myFun(LAMBDA, eigvals_F, clusterNum, F_eigvals_c1_pre)
%
% Inputs:
%   eigvals_F - c*1, c smallest eigenvalues of L_W
%   clusterNum - int, number of communities on the network
%   F_eigvals_sum_c_plus_pre - float, the sum of c+1 smallest eigenvalues of L_W (current value)
%
% Outputs:
%   update_w - bool, flag to indicate whether to update W
%   F_eigvals_sum_c_plus - float, the updated F_eigvals_sum_c_plus
%
% Author:  Y. Dong
% E-mail:  dyx_xjtu@163.com
% Created: 13 Feb, 2020

% Variables definition
ZR = 10e-11;
isWnotClear = true;
F_eigvals_sum_c = sum(eigvals_F(1:clusterNum));
F_eigvals_sum_c_plus = sum(eigvals_F(1:clusterNum+1));

% deterimate whether to update W
if F_eigvals_sum_c_plus_pre <= F_eigvals_sum_c_plus
    isWnotClear = false;
end

fprintf("f_c = %f; f_c_plus = %f \n", F_eigvals_sum_c, F_eigvals_sum_c_plus);