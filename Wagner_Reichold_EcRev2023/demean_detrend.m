function[ytilde,xtilde,x2tilde,x3tilde] = demean_detrend(y,x,type)
% demeaning or demeaning and detrending
%----------------------------------------------------------------------- 
%INPUTS:    y...        TxN matrix dependent variable
%           x...        TxN matrix integrated regressor
%           type...     demeaning (1) demeaning and detrending (2)  
%-----------------------------------------------------------------------
%OUTPUTS:   ytilde...   TxN-matrix transformed dependent variable
%           xtilde...   TxN-matrix transformed regressor
%           x2tilde...  TxN-matrix transformed squared regressor
%           x3tilde...  TxN-matrix transformed cubic regressor
%------------------------------------------------------------------------
% KR, April 2021
%------------------------------------------------------------------------
[T,~] = size(x);

if type == 1 % demeaning only
    D = ones(T,1);
elseif type == 2 % demeaning and detrending
    D = [ones(T,1),(1:T)'];
end

% projector on the orthogonal complement of D
P = eye(T) - D*((D'*D)\D');

% residuals of the projection
ytilde = P*y;
xtilde = P*x;
x2tilde = P*(x.^2);
x3tilde = P*(x.^3);

end



