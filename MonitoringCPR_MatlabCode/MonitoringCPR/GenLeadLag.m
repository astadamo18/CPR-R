%--------------------------------------------------------------------------
% function DeltaMat = GenLeadLag(y,n_lag,n_lead);
%
% Input:    y in R^{n \times T}
%           n_lag  ... max number of lags
%           n_lead ... max number of leads
%
% Output: DeltaMat in R^{n*(1+n_lag+n_lead) \times T}
%
%--------------------------------------------------------------------------
function DeltaMat = GenLeadLag(y,n_lag,n_lead);

[n,T1] = size(y);

% Initialization of results matrix:
lagmat = zeros(n*n_lag,T1);
leadmat = zeros(n*n_lead,T1);
DeltaMat = zeros(n*(1+n_lag+n_lead),T1);

% Creation of matrix with lags:
for j = 1:n_lag;
    lagmat((j-1)*n+1:j*n,j+1:end) = y(:,1:end-j);
end;

% Creation of matrix with leads:
for j = 1:n_lead;
    leadmat((j-1)*n+1:j*n,1:end-j) = y(:,j+1:end);
end;

DeltaMat = [y; lagmat; leadmat];

