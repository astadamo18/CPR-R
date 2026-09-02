%-------------------------------------------------------------------------
% function [Vmatrix] = Vhat_NEW(g)
%
% Function to compute the VCV Matrix of the coefficient estimates of 
% the augmented and partial summed regression (1) from Proc: IMOLS_NL
%
% Input:        g ... Regressor vector corresponding to (1)
%                     in R^{T \times #of variables in (1)
%
% Output: Vmatrix ... VCV of estimated of coefficients up to cond. LR-var.
%                   
% External procedures: none
%
% MW
%-------------------------------------------------------------------------
function [Vmatrix] = Vhat_NEW(g,A)

[T,size2] = size(g);

% Outer term:
%Term_1 = (1/T^0)*inv(A)*g'*g*inv(A);

% Central term:
% Summed quantities: 
S = cumsum(g);

SJ = zeros(T,size2);
SJ(2:end,:) = S(1:end-1,:);
ST2 = S(T,:);
ST = repmat(ST2,T,1);
DS = ST - SJ;

% % for j=1:T
% %     if j == 1;
% %         DS(j,:) = S(T,:);
% %     else;
% %     DS(j,:) = S(T,:)-S(j-1,:);
% %     end;
% % end;

%Term_2 = T^(-3)*inv(A)*DS'*DS*inv(A);

%Vmatrix = inv(Term_1)*Term_2*inv(Term_1);
gg = inv(g'*g);
%Vmatrix = gg\DS'*DS/gg;
Vh = gg*DS';
Vmatrix = Vh*Vh';