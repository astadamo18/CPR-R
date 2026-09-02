%----------------------------------------------------------
% function Z = gen_scmpr_regressors(x, powers)
%
% This function computes regressors for systems of cointegrating
% multivariate polynomial regressions, i.e., regressors Z such that:
%
%      y = cf_mat * Z + u,    Z = [z1, ..., zI]',    y = [y1, ..., yn]',
%
% with zi = t^i0 * x1^i1 * ... * xm^im.
%       
% Input:   x      ... an T x m matrix of integrated variables
%                     that do not do not cointegrate
%          powers ... an I x (m + 1) matrix containing the powers of
%                     the variables (t, x1, ..., xm); i-th row
%                     contains multi-index (i0, i1, ..., im)
%                      => i-th regressor: zi = t^i0 * x1^i1 * ... * xm^im
% Output:  Z      ... An T x I matrix containing the regressors for
%                     a system of cointegrating multivariate polynomial
%                     regressions, i.e., i-th column provides data for 
%                     regressor zi = t^i0 * x1^i1 * ... * xm^im
%
% External function:
%
% SV, April 2024
%----------------------------------------------------------
function Z = gen_scmpr_regressors(x, powers)

[T, m] = size(x);
[I, ~] = size(powers);

base = [(1:T)', x]; % base for exponents
Z = zeros(T, I);
for i = 1:I
    temp = 1;
    for j = 1:(m + 1)
        temp = temp .* base(:, j).^powers(i, j);
    end
    Z(:, i) = temp;
end

end
