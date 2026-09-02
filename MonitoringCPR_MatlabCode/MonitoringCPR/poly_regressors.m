%--------------------------------------------------------------------------
% function[Ft] = poly_regressors(T,p);
%
% Function to generate polynomial trend regressors
%
% Input  T ... Sample size
%        p ... maximal power of regressors: t^0,...,t^{p}
%
% Output Ft ... Regressors in R^{T \times (p+1)}
%
% TJV and MW, October 2009
%--------------------------------------------------------------------------
function[Ft] = poly_regressors(T,p);

if p == -1;    
    Ft = zeros(0,0);
elseif p > -1;
    Ft = zeros(T,p+1);
    one_vec = ones(T,1);
    time_vec = cumsum(one_vec);
    if p == 0;
        Ft = one_vec;
    elseif p == 1;
        Ft = [one_vec time_vec];
    else;
        Ft(:,1) = one_vec;
        Ft(:,2) = time_vec;
        for j=3:p+1;
            for s =1:T;
                Ft(s,j) = time_vec(s)^(j-1);
            end;
        end;
    end;
    
end;

