%------------------------------------------------------------------------------
%CTPUTests
% function result = FM_CPR(y,x,orders,w,deter,kern,band,deme,opt)
%
% Fully Modified OLS estimator for the cointegrating polynomial regression %(CPR) 
%(Procedure designed for the multivariate case; inference univariate right % now)
% 
% y = gamma*w + delta*deter + beta*X + u
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
%          w    ... Stationary variables in R^{T \times s}
%       deter   ... Deterministic variables, in R^{T \times q+1}
%                       (including constant)
%       kern    ... Specifies chosen kernel function:
%                tr ... Truncated
%                ba ... Bartlett
%                pa ... Parzen
%                bo ... Bohman
%                da ... Daniell
%                qs ... Quadratic Spectral
%       band    ... specifies the bandwidth chosen, 
%                   Either a number or Strings('And91','NW','AM92')
%       deme    ... Demeaning of residuals (1 yes, 0 no) in LR-variance
%                     computation (GENREICALLY WE TAKE 0 HERE, COULD BE
%                                   HARD WIRED INSIDE PROCEDURE)
%        opt (Use this when X is provided direcly)
%         .X      Matrix with all powers already in there: T times sum of vars.
%         .v      Delta(x)=v_t from the paper.
%         .Mxs    Correction terms without the Delta+ in front. (not A*)
%         .deter  Deterministic components
%         .P      e.g: X = [x1 x1^2 x2 x2^3], P = [0;2;4];
% Output:
%       result
%         .beta_fm      ... Coefficient matrix in R^{m \times 1}
%         .delta_fm     ... Coefficient matrix in R^{q \times 1}
%         .gamma_fm     ... Coefficient matrix in R^{s \times 1}
%         .t_beta       ... t-values for coeff.beta-FM (col vector)
%         .t_delta      ... t-values for coeff.theta-FM (col vector)
%         .t_gamma      ... t-values for coeff.gamma-FM (col vector)
%         .std_beta     ... std for coeff.beta-FM (col vector)
%         .std_delta    ... std for coeff.theta-FM (col vector)
%         .std_gamma    ... std for coeff.gamma-FM (col vector)
%         .u_plus       ... FM-OLS residuals in T^{m \times (T-1)}
%         .Omega_udotv  ... Cond. LR-var. based on OLS residuals
%         .varmat0      ... FM-OLS-HAC type VCV matrix for coeff. to
%                           stationary regressors (correct for inference!)
%         .varmat1      ... FM-OLS Variance Covariance matrix
%         .beta_ols     ... OLS coefficient matrix in R^{m \times 1}
%         .delta_ols    ... OLS coefficient matrix in R^{q \times 1}
%         .gamma_OLS    ... OLS coefficient matrix in R^{s \times 1}
%         .u_ols        ... OLS residuals in T^{m \times T}
%         .FInv         ... Untruncated (Z'Z)-inverse Matrix
%         .varmatOLS    ... OLS-HAC type ("wrong") variance covariance
%                           matrix for all coefficients (ignoring coint.)
%         .Fitted       ... Fitted Values
%
% External function: lr_varmod, lr_weights , lag, procedure that generates the included powers of the integrated regressors; procedure that % % % %       
% generates the necessary correction terms.
%
% Remark: COINT uses "deter"-detrending also for X_t!!!
%
% t-statistics for the stationary regressors have to be based on HAC covariance estimation
%----------------------------------------------------------
% (Required Files)
%    And_HAC91.m,  AndMon_HAC92.m,  AndMon_Stab.m, bwNW.m, GenCPRCorrVec.m,
%    GenPowerReg.m, GenVarPolyTerms.m, lag.m, lr_var.m, lr_weights.m,  var_m.m
%
%A) (2:end,:) egyszerű törlése ezeken a helyeken:
%u_ols(2:end,1) → u_ols (3 helyen fordul elő: And_HAC91, bwNW, lr_var/AndMon_HAC92 hívásokban)
%w(2:end,:) → w (2 helyen: az Astar(1:kw,:) sorban)
%y(2:end,:) → y (a yplus = ... sorban)
%Z(2:end,:) → Z (a bplus = iZZ*(...) sorban)

function result = FM_CPR(y,x,orders,w,deter,kern,band,deme,opt)
%% Basics

% We need a cell structure to accommodate all cases!!
%if nargin<9 || isempty(opt) 
    % In principle we could use correction factors from obs. nr. 2
    %ret = GenVarPolyTerms(x,orders);
    %X = ret.X;
    %P = ret.P;
    %Mstar = ret.Mstar; 

    %[T,m] = size(x);
    %v = diff(x);
    %javított rész
if nargin<9 || isempty(opt) 
    v = diff(x);            % ELŐSZÖR, a teljes (csonkolatlan) x-ből!
    y = y(2:end,:);          % harmonizálás: azonnali csonkolás
    x = x(2:end,:);
    if ~isempty(w)
        w = w(2:end,:);
    end
    deter = deter(2:end,:);

    ret = GenVarPolyTerms(x,orders);   % már a CSONKOLT x-en
    X = ret.X;
    P = ret.P;
    Mstar = ret.Mstar; 

    [T,m] = size(x);         % T mostantól már T0-1
else        % Providing Info Directly
    X = opt.X;
    P = opt.P;
    Mstar = opt.Mxs;
    deter = opt.deter;
    v = opt.v;
    
    m = length(P)-1;
    T = size(X,1);
end
kw = size(w,2); % kw is the number of additional stationary regressors
kd = size(deter,2);

J = [deter X];
Z = [w J];    % Regressors 
FInv = inv(Z'*Z); % Full (X'X)^-1 matrix!
 

%% OLS Regression
iZZ = inv(Z'*Z);
%b_ols = iZZ*Z_'*y(2:end);
b_ols = Z\y;
u_ols = y - Z*b_ols;
%% 
iww = inv(w'*w); % (w'w)^-1

%% FM Estimation (Hong & Wagner)
% (1) Constructing LR variance Estimators |||| LR means long-run
%if isequal(band,'And91')
%    bandw = And_HAC91([u_ols v],kern);
%elseif isequal(band,'NW')
%    bandw = bwNW([u_ols v],kern,0,[]);
%else
%    bandw = band;
%end
%if isequal(band,'AM92')
 %   [Lr,Dr,~] = AndMon_HAC92([u_ols v],kern,1,1,deme);
%else
%   [Lr,Dr,~] = lr_var([u_ols v],kern, bandw, deme);
%end
%disp('From FM_CPR'),Lr, Dr, Sr

%% FM Estimation (Hong & Wagner)
% (1) Constructing LR variance Estimators |||| LR means long-run
v_dm = v - mean(v);   % harmonizalas: v demeanelese a HAC/LR-variancia becslesbe (Wagner-Reichold 2023, Remark 5)
if isequal(band,'And91')
    bandw = And_HAC91([u_ols v_dm],kern);
elseif isequal(band,'NW')
    bandw = bwNW([u_ols v_dm],kern,0,[]);
else
    bandw = band;
end
if isequal(band,'AM92')
    [Lr,Dr,~] = AndMon_HAC92([u_ols v_dm],kern,1,1,deme);
else
    [Lr,Dr,~] = lr_var([u_ols v_dm],kern, bandw, deme);
end

Omega_udotv = Lr(1,1)-Lr(1,2:end)*inv(Lr(2:end,2:end))*Lr(2:end,1);
Lr_vvvu = inv(Lr(2:end,2:end))*Lr(2:end,1);
%Lambda0 = Dr(1,2:end)-Lr_vvvu'*Dr(2:end,2:end);
Lambda0 = Dr(2:end,1) - Dr(2:end,2:end)*Lr_vvvu; 

% (2) Constructing Correction Terms
for i = 1:m % m is the number of integrated regressors
    ind = P(i)+1:P(i+1);
    Mstar(ind) = Lambda0(i)*Mstar(ind); % This is the complete Mstar.
end
Astar = zeros(size(Z,2),1); % Astar is a Vector
if ~isempty(w) % This is the case, when stationary regressors are included
    Astar(1:kw,:) = (1/T)*w'*u_ols - (1/T)*w'*v*Lr_vvvu;
end
Astar(end-P(end)+1:end,:) = Mstar;  
result.Astar = Astar;
result.Lambda0 = Lambda0;
% (3) FM Estimator
yplus = y-v*Lr_vvvu;
bplus = iZZ*(Z'*yplus-Astar); % FM Estimator for CPR.

%% Outputs
result.beta_fm = bplus(kw+kd+1:end,:);
result.delta_fm = bplus(kw+1:kw+kd,:);
result.gamma_fm = bplus(1:kw,:);
result.beta_ols = b_ols(kw+kd+1:end,:);
result.delta_ols = b_ols(kw+1:kw+kd,:);
result.gamma_ols = b_ols(1:kw,:);

result.u_plus = yplus - Z*bplus;
result.y_plus = yplus;
% In the SUR code this is called: uhat_plus - since it is a residual, not
% an error term!!!

result.Fitted = Z*bplus;
% result.endo_corr = v*Lr_vvvu;
result.u_ols = u_ols;

% PROPER FM-OLS HAC TYPE INFERENCE FOR COEFFICIENTS TO STATIONARY
% REGRESSORS:
%S = Z_.*(result.u_plus*ones(1,size(Z,2)));
if ~isempty(w)
    
    S = w.*(result.u_plus*ones(1,size(w,2)));
    if isequal(band,'AM92')
        [SLr,~,~] = AndMon_HAC92(S,kern,1,1,deme);
    else
        [SLr,~,~] = lr_var(S,kern, bandw, deme);
    end
    %result.varmat0 = T*iww*SLr*iww; 
    result.varmat0 = T*iww*SLr*iww; 
    tvm = diag(result.varmat0).^.5; %.^. is used for the root.
    result.std_gamma = tvm(1:kw);
    result.t_gamma = result.gamma_fm./tvm(1:kw);
else
    result.std_gamma = [];
    result.t_gamma = [];
    result.varmat0 = [];
end
% To make results consistent with the SUR Computations: 
if isequal(band,'And91')
    bandb = And_HAC91(result.u_plus,kern);
elseif isequal(band,'NW')
    bandb = bwNW(result.u_plus,kern,0,[]);
else
    bandb = band;
end

if isequal(bandb,'AM92')
    [Omega_udotv1,~,~] = AndMon_HAC92(result.u_plus,kern,1,1,deme);
else
    [Omega_udotv1,~,~] = lr_var(result.u_plus,kern, bandb, deme);
end

result.varmat = Omega_udotv*inv(J'*J);
result.varmat1 = Omega_udotv1*inv(J'*J);
%result.varmat = Omega_udotv*inv(J(2:end,:)'*J(2:end,:));
%result.varmat1 = Omega_udotv1*inv(J(2:end,:)'*J(2:end,:)); 
%result.varmat1 = Omega_udotv1*inv(Z(2:end,:)'*Z(2:end,:)); 
tvm1 = diag(result.varmat).^.5; 
result.std_delta = tvm1(1:kd);
result.std_beta = tvm1(kd+1:end);
result.t_delta = result.delta_fm./tvm1(1:kd);
result.t_beta = result.beta_fm./tvm1(kd+1:end); 

result.Omega_udotv = Omega_udotv;
result.Omega_udotv1 = Omega_udotv1;
result.FInv = FInv; % FINV is just (Z'Z)^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFERENCE PART FOR OLS TYPE ANALYSIS (INCORRECT APPROACH - SINCE COINT. REGRESSION!!!):
% GIVES OLS HAC COVARIANCE MATRIX FOR ALL COEFFICIENTS
S_ols = Z.*(result.u_ols*ones(1,size(Z,2)));   
    if isequal(band,'AM92')
        [SLr_ols,~,~] = AndMon_HAC92(S_ols,kern,1,1,deme);
    else
        [SLr_ols,~,~] = lr_var(S_ols,kern,bandw,deme);
    end 
%result.varmatOLS = T*FInv*SLr_ols*FInv;
result.varmatOLS = T*FInv*SLr_ols*FInv;
result.SLr_ols = SLr_ols;

% If And91 is used the results from here differ from those in SUCPRnc,
% because there at this step the bandwidth is computed for "Z.*uhat_ols"
% Maybe change that here as well at some point.
return;

