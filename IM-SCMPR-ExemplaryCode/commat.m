%----------------------------------------------------------
% function result = commat(m, n)
%
% This function computes the commutation matrix for an m x n matrix, i.e.
% it computes the mn x mn matrix K such that:
%
%   K * vec(A) = vec(A')
%
% for any matrix A of dimension m x n. Here, vec denotes vectorization
% by columns.
%       
% Input:   m      ... number of rows of matrix that is vectorized
%          n      ... number of columns of matrix that is vectorized
% Output:  result ... the commuation matrix for an m x n matrix
%
% External functions:
%
% SV, April 2024
%----------------------------------------------------------

function result = commat(m, n)

mn = m * n;
temp = reshape(1:mn, m, n)';
result = eye(mn);
result = result(reshape(temp, mn, 1), :);

end
