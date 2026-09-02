%----------------------------------------------------------
% function [Omega,Delta,Sigma] = lr_varmod(u,kern,band,deme)
%
% This function computes the long-run variance Omega,
% the one sided long-run variance Delta (starting with Lag 0) and
% the variance Sigma from an input matrix of residuals.
% in a system of the form:
%       
% Input:   u     ... Residual matrix, in R^{T \times m}
%          kern  ... Specifies chosen kernel function:
%                tr ... Truncated
%                ba ... Bartlett
%                pa ... Parzen
%                bo ... Bohman
%                da ... Daniell
%                qs ... Quadratic Spectral
%          band  ... specifies the bandwidth chosen, integer
%                    in the range 1,...,T
%          deme  ... Demeaning of residuals (1 yes, 0 no)
% Output:  Omega ... Long-run variance matrix
%          Delta ... One-sided long-run variance matrix
%          Sigma ... Variance matrix
%
% External function: lr_weights;
%
% TV and MW, August 2009,
% modification by RK and OS, August 2014
% KR error fixed (line 79), April 2016.
%----------------------------------------------------------
function[Omega,Delta,Sigma] = lr_varmod(u,kern,band,deme)

[T,m] = size(u);

% Demeaning residuals (full vector demeaning):
if deme == 1
    u = u-ones(size(u,1),1)*mean(u);
end

% Call of weight vector generating function:
[w,j_max] = lr_weights(T,kern,band);

%w = ones(T-1,1);
%j_max = T-1;

% Initialization:
Omega = zeros(m,m); 
Delta = zeros(m,m);

% Variance (scaling only by 1/T): 
Sigma = 1/T*(u'*u);

if (isequal(kern,'qs') || isequal(kern,'da'))
    ws = [0;w];
    c = zeros(1,T);
    R = toeplitz(ws,c)';
    R(T,:) = zeros(1,T);
    %Delta = 1/T*(u'*R*u);
    Delta1 = R*u;
    Delta = 1/T*(u'*Delta1);
    Omega = Delta + Delta';    
else
    for j = 1:j_max
  
      % Residual summation (not using transposition equality):
      T1 = u(j+1:end,:)'*u(1:T-j,:)/T;
      T2 = u(1:T-j,:)'*u(j+1:end,:)/T;
    
      % Note: T1 and T2' sometimes minor numerical differences between 
      % T1 and T2' therefore summing both separately.
    
      % Scaling with weight function:
      T_Omega = w(j)*(T1+T2);
      T_Delta = w(j)*T1;
      % Summation:
      Omega = Omega + T_Omega;
      Delta = Delta + T_Delta;    
    
    end % j-Loop
    Delta = Delta';
    % This is the same as if we would define 
    % T_Delta = w(j)*T2; and Delta = Delta + T_Delta;
    % and delete line 79.
end

% Adding contemporaneous variance:
Omega = Omega + Sigma;
Delta = Delta + Sigma;
end