function result = DOLS_CPR(y,x,orders,deter,n_lead,n_lag,kern,band,deme,opt)
%------------------------------------------------------------------------------
% function result = DOLS_CPR(y,x,orders,deter,n_lead,n_lag,kern,band,deme,opt)
%
% Dynamic OLS estimator for the cointegrating polynomial regression (CPR) 
% (Procedure designed for the multivariate case; inference univariate right now)
% 
%      y = delta*deter + beta*X + chi*x_LeadLag + u
% 
% where X contains x and powers thereof and x_LeadLag contains leads and
% lags of first differences of x only (withoud powers)
%
% Input:   y    ... Dependent variable, in R^{T \times 1}
%          x    ... Explanatory variable, in R^{T \times m} 
%                   There needs to be at least one x variable!
%     orders    ... vectors that specify the orders of x_i included
%                   Does a cell array need to be defined? 
%                   - if [x^2 x^4] is needed, put in cellarray.
%                   (ex: orders{1}=[1 4]; orders{2}=[2 3 5]; etc)
%                   - if scalar, all x_i's have the same max orders
%                   - if vector, they are x_i's max orders
%       deter   ... Deterministic variables, in R^{T \times q+1}
%                       (including constant)
%        n_lead ... number of leading first differences for DOLS regression
%        n_lag  ... number of contemporaneous + lagging first differences for DOLS regression
%                   NOTE: if n_lag==0, then also no contemporaneous lag is
%                   included!
%                   if (n_lead + n_lag) == 0, then simply the original equation (augmented 
%                   by powers of the regresors) is estimated by OLS;
%         kern     ... Kernel for long-run variance estimation of resids.
%         band     ... Bandwidth choice for long-run variance estimation,
%                      number or 'AND91'; or 'NW', 'AM92'
%       deme    ... Demeaning of residuals (1 yes, 0 no) in LR-variance
%                     computation (GENERICALLY WE TAKE 0 HERE, COULD BE
%                                   HARD WIRED INSIDE PROCEDURE)
%        opt (Use this when X is provided direcly)
%         .X      Matrix with all powers already in there: T times sum of vars.
%         .v      Delta(x)=v_t from the paper.
%         .Mxs    Correction terms without the Delta+ in front. (not A*)
%         .deter  Deterministic components
%         .P      e.g: X = [x1 x1^2 x2 x2^3], P = [0;2;4];
% Output:
%       result
%         .beta_dols    ... Coefficient matrix in R^{(m + #powers of x) \times 1}
%         .delta_dols   ... Coefficient matrix in R^{q \times 1}
%         %.chi_dols     ... Coefficient matrix in R^{(1 + n_lead + n_lag) \times 1}
%                           ... format for leads and lags in regression [y; lagmat; leadmat]
%         .t_beta       ... t-values for coeff.beta_dols (col vector)
%         .t_delta      ... t-values for coeff.theta_dols (col vector)
%         %.t_chi        ... t-values for coeff.chi_dols (col vector)
%         .std_beta     ... std for coeff.beta_dols (col vector)
%         .std_delta    ... std for coeff.theta_dols (col vector)
%         %.std_chi      ... std for coeff.chi_dols (col vector)
%         .u_dols       ... DOLS residuals in R^{(T-(1+n_lead+n_lag) \times 1}
%         .Omega_udotv  ... method 1 (OLS): Cond. LR-var. based on OLS residuals
%                           method 2 (DOLS): LR-var. based on DOLS
%                           residuals
%         .varmat_dols  ... DOLS Variance Covariance matrix
%         .beta_ols     ... OLS coefficient matrix in R^{(m + #powers of x)  \times 1}
%         .delta_ols    ... OLS coefficient matrix in R^{q \times 1}
%         .u_ols        ... OLS residuals in R^{T \times 1}
%
% External function: lr_var,lr_weights,lag, procedure that generates the included 
% powers of the integrated regressors; procedure that generates the necessary correction terms.
%
% Remark: COINT uses "deter"-detrending also for X_t!!!
%
% t-statistics for the stationary regressors have to be based on HAC covariance estimation
%----------------------------------------------------------
% (Required Files)
%   GenPowerReg.m, GenCPRCorrVec.m, GenLeadLag.m, GenVarPolyTerms.mat, lr_varmod.m, 
%   lr_weights.m, AndMon_HAC92.m, And_HAC91.m, AndMon_Stab.m, var_m.m
%
% ###
% written (or better "modified", based upon MW's FM_CPR) by Stefan Schneeberger (s2),
% April 2011, IHS Vienna, stefan.schneeberger@gmx.at
% ###
% Modified by MW, April 2011

%% Basics

% We need a cell structure to accommodate all cases!!
if nargin<10 || isempty(opt) 
   ret = GenVarPolyTerms(x,orders);
   X = ret.X;
   P = ret.P;
   [T,m] = size(x);
   v = diff(x);
else        % Providing Info Directly
   X = opt.X;
   P = opt.P;
   deter = opt.deter;
   v = opt.v;

   m = length(P)-1;
   T = size(X,1);
end


kd = size(deter,2);
m_aug = size(X,2); %s2: number of regressors (including powers of x_is)
%size(deter)
%size(X)

Z = [deter X];    % Regressors
%size(Z)
%size(y)

%% OLS Regression
b_ols = Z\y;

u_ols = y - Z*b_ols;


%% DOLS Estimation (s2)

if ((n_lead+n_lag) == 0)
  % If neither leads nor lags are included: OLS estimation of original 
  % equation (computed above):
  
  Z_trunc = Z; % needed for varmat_dols computation later
  y_trunc = y;
  
  b_dols = b_ols;
  u_dols = u_ols;
  
else
  % (1) construct additional regressors (leads and lags of delta(x_i)):

  X_LeadLag = GenLeadLag(diff(x)',n_lag,n_lead)';

  % (2) construct (augmented and truncated) regressor matrix (including leads
  %     and lags of delta(x_i)) and truncate y (explanatory variable) accordingly:

  Z_untrunc = [Z(2:end,:) X_LeadLag];
  Z_trunc = Z_untrunc((n_lag+1):(end-n_lead),:);
  y_trunc = y((n_lag+2):(end-n_lead),:);

  % (3) run DOLS regression:
  
  b_dols = Z_trunc\y_trunc;
  u_dols = y_trunc - Z_trunc*b_dols;

end
%% Comment:
% Note that u_dols is computed based on the truncated versions of y and Z.
% That implies that (if plugged-in the lr_var-function via 'resid'), only
% the truncated values are used for the Long-Run Estimation
%% LR Estimation (s2)
% Constructing LR variance estimators

% version 1 (ols_residuals):
% version 2 (dols_residuals):
version = 2; % choose 1 or 2

if (version == 1)
  resid = [u_ols(2:end,1) v];  
else
  resid = u_dols;
end

if isequal(band,'And91');
   bandw = And_HAC91(resid,kern);
    [Lr,Dr,Sr] = lr_varmod(resid,kern, bandw, deme);
elseif isequal(band,'NW');
   bandw = bwNW(resid,kern,0,[]);
    [Lr,Dr,Sr] = lr_varmod(resid,kern, bandw, deme);
elseif isequal(band,'NWfake');
   bandw = floor(4*(length(y_trunc)/100)^(2/9));
   [Lr,Dr,Sr] = lr_varmod(resid,kern,bandw,deme);
elseif isequal(band,'AM92');
    [Lr,Dr,Sr] = AndMon_HAC92(resid,kern,1,1,deme);
    bandw = band;
else
   bandw = band;
   [Lr,Dr,Sr] = lr_var(resid,kern, bandw, deme);
end

% if isequal(band,'AM92');
%    [Lr,Dr,Sr] = AndMon_HAC92(resid,kern,1,1,deme);
% else
%    [Lr,Dr,Sr] = lr_varmod(resid,kern, bandw, deme);
% end

if (version == 1)
   % compute conditional LR variance estimator
   Omega_udotv = Lr(1,1)-Lr(1,2:end)*inv(Lr(2:end,2:end))*Lr(2:end,1);
else
   Omega_udotv = Lr;
end


%% Outputs (s2)
result.beta_dols = b_dols(kd+1:kd+m_aug,:);
result.delta_dols = b_dols(1:kd,:); 
result.chi_dols = b_dols(kd+m_aug+1:end,:); 
% chi (= coeff of add. leads and lags) not interesting, hence
% results skipped from output;

result.beta_ols = b_ols(kd+1:end,:);
result.delta_ols = b_ols(1:kd,:);

result.Fitted = Z*b_dols(1:(kd+m_aug));

result.u_ols = u_ols;
result.u_dols = u_dols;


result.varmat_dols = Omega_udotv * inv(Z_trunc'*Z_trunc);
tvm_dols = diag(result.varmat_dols).^.5;

result.std_delta = tvm_dols(1:kd);
result.std_beta = tvm_dols((kd+1):(kd+m_aug));
%result.std_chi = tvm_dols((kd+m_aug+1):end); not interesting

result.t_delta = result.delta_dols./tvm_dols(1:kd);
result.t_beta = result.beta_dols./tvm_dols((kd+1):(kd+m_aug));
%result.t_chi = result.chi_dols./tvm_dols((kd+m_aug+1):end); not interesting

result.Omega_udotv = Omega_udotv;

return;