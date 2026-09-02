%--------------------------------------------------------------------------
%CTPUTests
% function corr_term = GenCPRCorrVec(x,all,powvec);
%
% Function to return a vector of correction terms corresponding to power of   
% (scalar) integrated regressor: Generate CPR Correction Vector.
% [CPR = Cointegrating Polynomial Regression] 
%
% Input:            x ... (Single) Regressor in R^{T \times 1}
%                 all ... If all == yes, then all powers from 1 to powvec 
%                         - which then is scalar - are returned.
%                         If all == no, then the powers given in powvec
%                         - then a vector - are returned.
%              powvec ... Vector in R^{1 \times NumPow}
%
% Output:   corr_term ... Matrix in R^{NumPow \times 1} with corr. terms 
%                         corresponding to powers (1: T, 2: 2\sum x_t,...)
%--------------------------------------------------------------------------
function corr_term = GenCPRCorrVec(x,all,powvec);

% In the computation all power terms from 1 to highest power need to be
% computed anyway. 

[T,sizex2] = size(x);
if isempty(powvec)
    corr_term = [];
    return
end
powvec = powvec(:)';    % (2010.11.8) Added by Hong

if isequal(all,'yes');
    numpow = powvec;
    maxpower = powvec;
elseif isequal(all,'no');
    numpow = size(powvec,2);
    maxpower = powvec(numpow);
end;

% Computation of powers of regressors:
sum_matrix = ones(T,numpow);

sum_matrix(:,2) = x;
for j = 3:maxpower;
    sum_matrix(:,j) = sum_matrix(:,j-1).*x;
end;
    
% Computation of sums of powers:
for j = 1:maxpower
    sum_vec(j,1) = sum(sum_matrix(:,j));
end;
    
% Multiplication with "p-1" scalar:
prevec = ones(maxpower,1);
prevec = cumsum(prevec);

full_vec = prevec.*sum_vec;

% Selection of relevant elements:
if isequal(all,'yes');
    corr_term = full_vec;
elseif isequal(all,'no');
    for j = 1:numpow;
        corr_term(j,1) = full_vec(powvec(j),1);
    end;
end;

