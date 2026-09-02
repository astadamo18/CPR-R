%-------------------------------------------------------------------------
% function [theta_hat1,u_hat1,Vhat,theta_hat2,u_hat2]=IMOLS_NL(y,deterc,xL,xNL,augtype,selector,Vcalc)
% 
% Function to compute the "first and second step" OLS regression of the IM-OLS
% estimator for nonlinear, linear in parameters cointegrating relationships: 
%               y = D'*delta + xL'*betaL + xNL'*betaNL + u (*)
%  Step 1 Partial Sum - with linear augmentation: 
%           Sy = SD'*delta + SxL'*betaL + SxNL'*betaNL + 
%                + xL'*gamma + Su       (1)      
%  Step 1 Partial Sum - with full augmentation:
%           Sy = SD'*delta + SxL'*betaL + SxNL'*betaNL +
%                +xL'*gamma_1 + xNL*gamma2 + Su            (1)
%  Step 2: On either of the two versions above - augmentation by Z 
%          (using all nondet. regressors): x* corresponds to augmentation
%          type (= xL or (xL,xNL))
%
%     Sy = SD'*delta + SXall'*beta + x*'*gamma + Z'*lambda + u    (2)
% -------------------------------------------------------------------------
% Input     y ... Dependent variable in R^T
%           D ... Deterministic components in R^{T x dim(D)}
%          xL ... Linear Integrated regresors in R^{T \times dim(xL)}
%         xNL ... Non-linear Transformations of Integrated regresors in
%                 R^{T \times dim(xNL)}, i.e. in polynomial case: [x_t^2,...,x_t^p]
%     augtype ... Augmentation type
%                 1 = Linear Augmentation
%                 2 = Full Augmentation
%    selector ... 1 = Only regression (1)
%                 2 = Only regression (2)
%                 3 = Both regressions (1) and (2)
%       Vcalc ... 0 = Vhat not computed
%                 1 = Vhat is computed
%
% Output:   theta_hat1 ... [delta',beta',gamma']' coefficient vector reg (1)
%               u_hat1 ... Residual vector in R^T of reg (1)
%                 Vhat ... Estimated VCV (up to Lambda^2) of reg (1)
%           theta_hat2 ... [delta',beta',gamma',lambda']' coefficient vector reg (2)
%               u_hat2 ... Residual vector in R^T of reg (2)
%
% External functions: poly_regressors, Z_regressors, Vhat_Reg1
%------------------------------------------------------------------------- 
function result =IMOLS_NL(y,deterc,xL,xNL,augtype,selector,Vcalc)


% Sample sizes and variable numbers: 
[T,y1] = size(y);
[TxL,kL] = size(xL);
[TxNL,kNL] = size(xNL);
[TxD, kD] = size(deterc);
xALL = [xL xNL];        % All stochastic regressors

% Generate Partial Summed Regressors:
  Sy = cumsum(y);
  SxALL = cumsum(xALL);
  SD = cumsum(deterc);

% Generation of regressor matrix corr. to (1):
  if augtype == 1;
      Xmat = [SD SxALL xL];
      %Xstochmat = [SxALL xL];
  elseif augtype == 2;
      Xmat = [SD SxALL xALL];
      %Xstochmat = [SxALL xALL];
  end;
  
  
if selector ~= 2;    
    % Regression (1) itself:    
        % Coefficient vector:
        theta_hat1 = Xmat\Sy;
        result.theta_hat1 = theta_hat1;
        % Residuals (in case needed):
        result.u_hat1 = Sy - Xmat*theta_hat1;
        
        %Fitted Values:
        result.Fitted = [deterc xL xNL]*theta_hat1(1:(kD+kL+kNL),1);
    
    if Vcalc == 0;
        result.Vhat = [];
    elseif Vcalc == 1;
        % Computation of Vmatrix: 
        result.Vhat = Vhat_NEW(Xmat);
    end;
    
end;

if selector ~= 1;

    % Regression (2) itself:

    % Generation of Z-regressors (based on summed regressors):
    % Automatically corresponding to augmentation style - since all regressors
    % in (1) will be Z-double summed.
    Z_reg = Z_regressors(Xmat);
    X2mat = [Xmat Z_reg];
        % Coefficient vector:
        theta_hat2 = X2mat\Sy;
        result.theta_hat2=theta_hat2;
        % Residual vector: 
        result.u_hat2 = Sy - X2mat*theta_hat2;
  
end;

if selector == 1;
    result.theta_hat2 = [];
    result.u_hat2 = [];
elseif selector == 2;
    result.theta_hat1 = [];
    result.u_hat1 = [];
    result.Vhat = [];
end;
