%----------------------------------------------------------
% function [coeff, fitted, resid, coeff_var, resid_mod, coeff_var_mod, b, coeff_r, coeff_f, coeff_var_f] = im_scmpr(y, x, powers, D, d, weights, kern, band)
%
% This function computes IM-OLS and IM-GLS estimates for a system of 
% cointegrating multivariate polynomial regressions. Consider:
%
%      y = cf_mat * Z + u,    Z = [z1, ..., zI]',    y = [y1, ..., yn]',
%
% with zi = t^i0 * x1^i1 * ... * xm^im. IM regression:
%
%      Sy = cf_mat * SZ + cf_aug * x + Su,    x = [x1, ..., xm]',
%
% where vec([cf_mat, cf_aug]') = D * cf_free + d.
%       
% Input:   y             ... an T x n matrix of integrated variables 
%                            that cointegrate with x; left-hand side
%                            of the equation system; dependent variables
%          x             ... an T x m matrix of integrated variables
%                            that do not do not cointegrate
%          powers        ... an I x (m + 1) matrix containing the powers of
%                            the variables (t, x1, ..., xm); i-th row
%                            contains multi-index (i0, i1, ..., im) such
%                            that i-th regressor is given by:
%                            zi = t^i0 * x1^i1 * ... * xm^im
%          D             ... an (I + m)n x g restriction matrix such that:
%                            vec([cf_mat, cf_aug]') = D * cf_free + d
%                            Default is (I + m)n x (I + m)n identity matrix
%          d             ... an (I + m)n x 1 restriction vector (see D)
%                            Default is (I + m)n x 1 vector of zeros
%          weights       ... the weightning matrix for IM-GLS. Either:
%                            an n x n matrix ... matrix to use
%                            'PaOg1'         ... uses hat(Omega)_uu^(-1)
%                            'PaOg2'         ... uses hat(Omega)_u.v^(-1)
%                            Default is n x n identity matrix
%          kern          ... specifies chosen kernel function:
%                            'tr' ... Truncated
%                            'ba' ... Bartlett
%                            'pa' ... Parzen
%                            'bo' ... Bohman
%                            'da' ... Daniell
%                            'qs' ... Quadratic Spectral
%                            Default is 'ba'
%          band          ... specifies bandwith chosen. Either:
%                            int in {1, ..., T} ... bandwidth to use 
%                            'And'              ... data dependent rule,
%                                                   Andrews (1991)
%                            'NW'               ... data dependent rule,
%                                                   Newey & West (1987)
%                            'NWT'              ... sample size dependent 
%                                                   rule of thumb
%                            Default is 'And'
% Output:  coeff         ... an n x (I + m) matrix of unrestricted IM-OLS
%                            estimates of [cf_mat, cf_aug]
%          fitted        ... an T x n matrix of fitted values from 
%                            unrestricted IM-OLS estimation
%          resid         ... an T x n matrix of residuals from 
%                            unrestricted IM-OLS estimation
%          coeff_var     ... an (I + m)n x (I + m)n matrix giving the 
%                            estimated 'asymptotic' covariance matrix of
%                            of vec(coeff') for standard-asymptotic theory
%          resid_mod     ... an T x n matrix of modified IM-OLS residuals
%                            used for fixed-b long-run covariance matrix
%                            estimation
%          coeff_var_mod ... an (I + m)n x (I + m)n matrix giving the 
%                            estimated 'asymptotic' covariance matrix of
%                            of vec(coeff') for fixed-b asymptotic theory
%          b             ... bandwidth-to-sample-size ratio computed using
%                            band and, if band is not int, diff(resid_mod).
%                            If band is int, then b = band/T.
%          coeff_r       ... an (I + m)n x 1 vector; the restricted IM-GLS
%                            estimates of vec([cf_mat, cf_aug]') with
%                            vec([cf_mat, cf_aug]') = D * cf_free + d
%          coeff_f       ... an g x 1 vector; the IM coefficient estimates
%                            cf_free
%          coeff_var_f   ... an g x g matrix giving the estimated 
%                            'asymptotic' covariance matrix of 
%                            coeff_f for standard-asymptotic theory
%
% External functions: lr_weights; lr_var; bwNW; And_HAC91;
%
% REMARKS: This function does not check whether the asymptotic rank
%          conditions hold, or, whether full design applies. In addition,
%          it does not check, whether D and d are correctly specified,
%          i.e., that the restrictions do not involve elements of cf_aug. 
%
% SV, April 2024
%----------------------------------------------------------
function [coeff, fitted, resid, coeff_var, resid_mod, coeff_var_mod, b, coeff_r, coeff_f, coeff_var_f] = im_scmpr(y, x, powers, D, d, weights, kern, band)

[T, n] = size(y);
[~, m] = size(x);
[I, ~] = size(powers);

if nargin < 8
    band = 'And';
    if nargin < 7
        kern = 'ba';
        if nargin < 6
            weights = eye(n);
            if nargin < 5
                d = zeros((I + m) * n, 1);
                if nargin < 4
                    D = eye((I + m) * n);
                end
            end
        end
    end
end

base = [(1:T)', x]; % base for exponents
Z = zeros(T, I);

for i = 1:I
    temp = 1;
    for j = 1:(m + 1)
        temp = temp .* base(:, j).^powers(i, j);
    end
    Z(:, i) = temp;
end

% 1) Unrestricted case:

SZ = [cumsum(Z), x];
Sy = cumsum(y);

coeff = SZ\Sy;
fitted = SZ * coeff;
resid = Sy - fitted;
coeff = coeff';

S2Z = cumsum(SZ);
C = S2Z(end, :) - S2Z;
C = [S2Z(end, :)', C(1:end-1, :)']';

SZ_cross = SZ' * SZ;
SZ_cross_inv = inv(SZ_cross);
calM_hat = SZ_cross_inv * (C' * C) * SZ_cross_inv;

% Conventional asymptotic theory:

coeff_ols = Z\y;
resid_ols = y - Z * coeff_ols;
resid_convl = [resid_ols(2:end, :), diff(x)];

if isequal(band, 'And')
    bandw = And_HAC91(resid_convl, kern);
elseif isequal(band, 'NW')
    bandw = bwNW(resid_convl, kern, 0, []);
elseif isequal(band, 'NWT')
    bandw = floor(4 * (T/100)^(2/9));
else
    bandw = band;
end

[lrvar, ~, ~] = lr_var(resid_convl, kern, bandw, 0);

yidx = 1:n;
xidx = n + (1:m);

temp = inv(lrvar(xidx, xidx)) * lrvar(xidx, yidx);
lrvar_cond = lrvar(yidx, yidx) - lrvar(yidx, xidx) * temp;

coeff_var = kron(lrvar_cond, calM_hat);

% Fixed-b asymptotic theory:

M = cumsum(C);
SZ_mod = [SZ, M];
resid_mod = Sy - SZ_mod * (SZ_mod\Sy);
d_resid_mod = diff(resid_mod);

if isequal(band, 'And')
    bandw = And_HAC91(d_resid_mod, kern);
elseif isequal(band, 'NW')
    bandw = bwNW(d_resid_mod, kern, 0, []);
elseif isequal(band, 'NWT')
    bandw = floor(4 * (T/100)^(2/9));
else
    bandw = band;
end

b = bandw/T;

[lrvar_cond_mod, ~, ~] = lr_var(d_resid_mod, kern, bandw, 0);

coeff_var_mod = kron(lrvar_cond_mod, calM_hat);

% 2) Restricted case:

if isequal(weights, 'PaOg1')
    W_hat = inv(lrvar(yidx, yidx));
elseif isequal(weights, 'PaOg2')
    W_hat = inv(lrvar_cond);
else
    W_hat = weights;
end

A_hat = D' * kron(W_hat, SZ_cross) * D;
B_hat = D' * kron(W_hat * lrvar_cond * W_hat, C' * C) * D;
temp = SZ' * Sy * W_hat;
temp = D' * (temp(:) - kron(W_hat, SZ_cross) * d);

coeff_f = A_hat\temp;
coeff_r = D * coeff_f + d;

A_hat_inv = inv(A_hat);

coeff_var_f = A_hat_inv * B_hat * A_hat_inv;

end
