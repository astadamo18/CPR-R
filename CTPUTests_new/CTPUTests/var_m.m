%-----------------------------------------------------------
%CTPUTests
% function [coeffs,resids,aic,bic] = var_m(y,p)
%
% Function to estimate VAR with OLS, truncation
% of observations lost by lagging!
% Normalization: a_0 = I,
% y_t = a1*y_t-1 + ... + a_p + e_t
%
% Input:         y ... Observations, in R^{T \times s}
%                p ... Order of the VAR
% Output:   coeffs ... [a_1,...,a_p]
%           resids ... Residuals in R^{T-p \times s}
%              aic ... Akaike Information Crit. A
%              bic ... Akaike Information Crit. B
% External Functions: lag, trimr
% NOTE: No deterministic variables included in this version!
%-----------------------------------------------------------
function [coeffs,resids,aic,bic] = var_m(y,p)

[T,s] = size(y);

%Generate lagged regressor matrix
for i = 1:p,
    regs(:,(1+(i-1)*s):i*s) = lag(y,i);    
end;
%regs = mlag(y,p);
%Cut lost observations
y_eff = trimr(y,p,0);
r_eff = trimr(regs,p,0);
%Regression
coeffs = inv(r_eff'*r_eff)*(r_eff'*y_eff);
resids = y_eff - r_eff*coeffs;
%Order estimation AIC and BIC
VCV = resids'*resids;
rVCV = VCV/(T-p*s);

%AIC
aic = log(det(rVCV)) + (2*p*s*s)/T;
%BIC 
bic = log(det(rVCV)) + (p*s*s*log(T))/T;
%Coefficients as estimated above are "transposed", thus:
coeffs = coeffs';


