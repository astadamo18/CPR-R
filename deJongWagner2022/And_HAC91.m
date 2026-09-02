%--------------------------------------------------------------------------
% function [bandwidth] = And_HAC91(v,kern)
%
% Function to implement the automatic bandwidth selection of
% Andrews (1991). Here the AR(1) individual version is implemented.
%
% Input:         v ... Data in R^{T \times dim(v)}
%             kern ... Kernel function: 
%                      'tr' ... Truncated
%                      'ba' ... Bartlett
%                      'pa' ... Parzen
%                      'th' ... Tukey-Hanning
%                      'qs' ... Quadratic-Spectral
%
% Output:  bandwidth ... Bandwidth selected.
%
% REMARK: Tukey-Hanning kernel not implemented elsewhere as of now
% (25/12/2009)
%--------------------------------------------------------------------------
function [bandwidth] = And_HAC91(v,kern);

[T,dimv] = size(v);

rhovec = zeros(dimv,1);
sigma2vec = zeros(dimv,1);

% Computation of rho and sigma^2 for each coordinate of v:
for j = 1:dimv;
    y = v(2:end,j);
    x = v(1:(end-1),j);
    rhovec(j,1) = x\y;
    sigma2vec(j,1) = (1/T)*(y-x*rhovec(j,1))'*(y-x*rhovec(j,1));
end;

% Computation of alpha(2): And91, 835, (6.4)
denom = 0;
for k = 1:dimv;
    denom = denom + (sigma2vec(k,1)^2)/((1-rhovec(k,1))^4);
end;
numer2 = 0;
for k = 1:dimv;
    numer2 = numer2 + (4*rhovec(k,1)^2*sigma2vec(k,1)^2)/((1-rhovec(k,1))^8); 
end;
a2 = numer2/denom;



% Computation of alpha(2): And91, 835, (6.4)
numer1 = 0;
for k = 1:dimv;
    numer1 = numer1 + (4*rhovec(k,1)^2*sigma2vec(k,1)^2)/((1-rhovec(k,1))^6*(1+rhovec(k,1))^2); 
end;
a1 = numer1/denom;

% Computation of bandwidth: And91, 834, (6.2) and footnote 5.
if kern == 'tr';
    bandwidth = 0.661*(a2*T)^(1/5);
elseif kern == 'ba';
    bandwidth = 1.1447*(a1*T)^(1/3);
elseif kern == 'pa';
    bandwidth = 2.6614*(a2*T)^(1/5);
elseif kern == 'th';
    bandwidth = 1.7462*(a2*T)^(1/5);
elseif kern == 'qs';
    bandwidth = 1.3221*(a2*T)^(1/5);
end;

bandwidth = ceil(bandwidth);
end