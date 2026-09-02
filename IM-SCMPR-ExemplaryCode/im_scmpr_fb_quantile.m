%----------------------------------------------------------
% function result = im_scmpr_fb_quantile(n, powers, hypo, b, kern, p, T, repetitions)
%
% This function computes (approximate) quantiles of the fixed-b limiting
% distribution of mod. IM-OLS residuals based Wald- and t-type statistics
% for systems of cointegrating multivariate polynomial regressions
%       
% Input:   n           ... number of equations, a positive integer
%          powers      ... an I x (m + 1) matrix containing the powers of
%                          the variables (t, x1, ..., xm); i-th row
%                          contains multi-index (i0, i1, ..., im)
%                          => i-th regressor is t^i0 * x1^i1 * ... * xm^im
%          hypo        ... number of hypotheses per equation
%          b           ... bandwidth-to-sample-size ratio, with 0 < b <= 1
%          kern        ... Specifies chosen kernel function:
%                      'tr' ... Truncated
%                      'ba' ... Bartlett
%                      'pa' ... Parzen
%                      'bo' ... Bohman
%                      'da' ... Daniell
%                      'qs' ... Quadratic Spectral
%          p           ... the p-quantile(s) to compute, 0 < p <= 1
%          T           ... sample size for approximating limiting
%                          quantities
%                          Default is 1000
%          repetitions ... number of repetitions for simulation
%                          Default is 10000
% Output:  wald        ... p-quantile(s) of fixed-b Wald-type statistic
%          t           ... p-quantile(s) of fixed-b t-type statistic
%                          Only available if hypo = n = 1, otherwise NaN
%
% External function: lr_weights; lr_var;
%
% SV, April 2024
%----------------------------------------------------------
function [wald, t] = im_scmpr_fb_quantile(n, powers, hypo, b, kern, p, T, repetitions)

if nargin < 8
    repetitions = 10000;
    if nargin < 7
        T = 1000;
    end
end

[I, m] = size(powers);
m = m - 1;

stats = zeros(repetitions, 2);

for k = 1:repetitions
    y = normrnd(0, 1, [T, n]);
    x = cumsum(normrnd(0, 1, [T, m]));
    base = [(1:T)', x]; % base for exponents
    Z = zeros(T, I);
    for i = 1:I
        temp = 1;
        for j = 1:(m + 1)
            temp = temp .* base(:, j).^powers(i, j);
        end
        Z(:, i) = temp;
    end
    SZ = [cumsum(Z), x];
    M = cumsum(SZ(T:-1:1, :));
    M = cumsum(M(T:-1:1, :));
    SZ_mod = [SZ, M];
    Sy = cumsum(y);
    resid_mod = Sy - SZ_mod * (SZ_mod\Sy);
    [Q, ~, ~] = lr_var(diff(resid_mod), kern, b * T, 0);
    calZ = normrnd(0, 1, [hypo * n, 1]);
    cent = kron(inv(Q), eye(hypo));
    stats(k, 1) = calZ' * cent * calZ;
    if hypo == 1 && n == 1
        stats(k, 2) = calZ * sqrt(cent);
    else
        stats(k, 2) = NaN;
    end
end

stats = sort(stats);
result = stats(ceil(repetitions * p), :);
wald = result(:, 1);
t = result(:, 2);

end
