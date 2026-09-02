%--------------------------------------------------------------------------
% function [Omega,Delta,Sigma] = AndMon_HAC92(u,kern,pw_lag,stab,deme)
%
% Function to compute VAR pre-whitened long-run variance estimates
% Combination of VAR pre-whitening with Andrews (1991) long-run variance
% estimation
%
% Input:   u     ... Residual matrix, in R^{T \times m}
%          kern  ... Specifies chosen kernel function:
%                      'tr' ... Truncated
%                      'ba' ... Bartlett
%                      'pa' ... Parzen
%                      'th' ... Tukey-Hanning
%                      'qs' ... Quadratic-Spectral
%         pw_lag ... Pre-whitening lag (1,2,3,...)
%           stab ... Stabilization of eigenvalues of a(1)^(-1) at 0.97.
%                    (0 = no; 1 = yes);
%           deme ... Demeaning of residuals (1 yes, 0 no)
% Output:  Omega ... Long-run variance matrix
%          Delta ... One-sided long-run variance matrix
%          Sigma ... Variance matrix
%
% External function: lr_weights; lr_var; And_HAC91, AndMon_Stab, var_m
% To be decided: correction factors or not?
%--------------------------------------------------------------------------
function[Omega,Delta,Sigma] = AndMon_HAC92(u,kern,pw_lag,stab,deme)

[T,m] = size(u);

% Demeaning residuals (full vector demeaning):
if deme == 1;
    u = u-ones(size(u,1),1)*mean(u);
end;

% VAR-prewhitening:
[coeffs,resids,aic,bic] = var_m(u,pw_lag);

% Compute coefficient polynomial evaluated at 1:
if pw_lag == 1;
    a1 = eye(m)-coeffs;
elseif pw_lag > 1;
    coeff_ext = [eye(m) coeffs];
    mult1 = kron(-ones(pw_lag,1),eye(m));
    multmat = [eye(m); mult1];
    a1 = coeff_ext*multmat;
end;

inva1 =inv(a1);

% Eigenvalue stabilization:
if stab == 1;
    inva1 = AndMon_Stab(inva1);
end;

% Andrews (1991) lag length computation:
band_pw = And_HAC91(resids,kern);

% Computation of variance: (Computed directly from residuals; no link to
% other parts of computations)
Sigma = 1/T*(u'*u);

% Long-run variance computation for pre-whitened residuals:
[Omega_pw,Delta_pw,Sigma_pw] = lr_var(resids,kern,band_pw,deme);

% "Re-coloring" of long-run variance:
Omega = inva1*Omega_pw*inva1';

% Computation of one-sided long-run covariance matrix:
% OPEN ISSUE / OPEN ISSUE
Lambda_pw  = Delta_pw  - Sigma_pw;
Delta = Sigma+inva1*Lambda_pw*inva1'+inva1*coeffs*Sigma;

