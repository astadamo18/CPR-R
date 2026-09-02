%--------------------------------------------------------------------------
% function[AIC_lag,AIC_lead,BIC_lag,BIC_lead] =
% LagLeadOrdersDOLS_CPR(y,x,degree,deter,max_lag,max_lead,symmet)
%
% Function to compute the optimal Lag and Lead Orders for DOLS as
% discussed in Choi and Kurozumi (2012).
%
% Input:     y     ... Dependent variable, in R^{n \times T}
%            x     ... Explanatory variable, in R^{m \times T}
%        degree    ... Degree of the last regressor (only the last regressor
%                      has higher order powers).
%         deter    ... Deterministic variables, in R^{q \times T}
%       max_lag    ... Max. number of Lags included.
%      max_lead    ... Max. number of Leads included.
%        symmet    ... if set == 1: then only same number of leads and lags
%                      is considered and searched over.
%                      Max. searched over is min(max_lag,max_lead)
%
% Output:   AIC_lag ... Opt. number of lags with AIC
%          AIC_lead ... Opt. number of leads with AIC
%           BIC_lag ... Opt. number of lags with BIC
%          BIC_lead ... Opt. number of leads with BIC
%
% External function: GenLeadLag
%
% MW, May 2010.
% 
% ###
% Modified "LagLeadOrdersDOLS.m" for Monitoring CPR (special case where only 
% the last regressor has higher order powers), by FK, July 2019.
%--------------------------------------------------------------------------
function[AIC_lag,AIC_lead,BIC_lag,BIC_lead] = ...
    LagLeadOrdersDOLS_CPR(y,x,degree,deter,max_lag,max_lead,symmet)

[n,T] = size(y);
[m,T] = size(x);
[q,T] = size(deter);

% Last regressor
x_k = x(m,:);

% Create regressor matrix:
if degree == 1
    Z = [deter; x];
else
    Z = [deter; x; repmat(x_k, [degree - 1, 1]).^(repmat((2:degree)', [1, T]))];
end

% Generate Delta x (and cut first observation):
DeltaX = x - lag(x',1)';
DeltaX = DeltaX(:,2:end);

% Cut first obs (lost due to first diff): 
Zs = Z(:,2:end);
ys = y(:,2:end);


% Computation without symmetry of lag and lead order imposed:
if symmet == 0;

%Initialization of results table:
AIC_table = zeros(max_lag+1,max_lead+1);
BIC_table = zeros(max_lag+1,max_lead+1);

for lag_num = 0:max_lag; % Loop over Lags:
    for lead_num = 0:max_lead; % Loop over Leads:
        
        % If no leads and lags are included: OLS on all data:
        if (lag_num + lead_num) == 0;
        OLSall = (Z'\y')';
        u_ols = y - OLSall*Z;
        act_samp = T;
        SSR = u_ols*u_ols';
        else;
        % DOLS estimation required:              
        % Generation of Leads and Lagged differences of regressors:
        dx_all = GenLeadLag(DeltaX,lag_num,lead_num);
        % Appropriate truncation of all regressors:
        AllUntrunc = [Zs; dx_all];
        AllTrunc = AllUntrunc(:,lag_num+1:end-lead_num);
        ytrunc = ys(:,lag_num+1:end-lead_num);
        % D-OLS estimation itself:
        DOLSall = (AllTrunc'\ytrunc')';
        % D-OLS residuals:
        u_dols = ytrunc - DOLSall*AllTrunc;
        act_samp = T-1-lag_num-lead_num;
        SSR = u_dols*u_dols';
        end;
        
        % Computation of information criteria:
        AIC_table(lag_num+1,lead_num+1) = act_samp*log(SSR/act_samp)+2*(m*(lag_num+lead_num+2)+2);
        BIC_table(lag_num+1,lead_num+1) = act_samp*log(SSR/act_samp)+log(act_samp)*(m*(lag_num+lead_num+2)+2);
        
    end; % Loop over Leads
end; % Loop over Lags

% Computation of optimal lead and lag orders:

    % AIC:
    A_long = AIC_table(:);
    [Asort,Aindex] = sort(A_long,'ascend');
    AIC_lead = ceil(Aindex(1,1)/(max_lag+1))-1;
    AIC_lag = mod(Aindex(1,1),(max_lag+1))-1;
    if AIC_lag == -1;
        AIC_lag = max_lag;
    end;
    % BIC:
    B_long = BIC_table(:);
    [Bsort,Bindex] = sort(B_long,'ascend');
    BIC_lead = ceil(Bindex(1,1)/(max_lag+1))-1;
    BIC_lag = mod(Bindex(1,1),(max_lag+1))-1;
    if BIC_lag == -1;
        BIC_lag = max_lag;
    end;

   
    %............Symmetric Model..................
% Computation for symmetric lag and lead orders:
elseif symmet == 1;
    
% Max order is minimum of prescribed values:
max_lag = min(max_lag,max_lead);
    
%Initialization of results table:
AIC_table = zeros(max_lag+1,1);
BIC_table = zeros(max_lag+1,1);

    
    for lag_comp = 0:max_lag;
        
        if lag_comp == 0;
        OLSall = (Z'\y')';
        u_ols = y - OLSall*Z;
        act_samp = T;
        SSR = u_ols*u_ols';
        else;
        % DOLS estimation required:              
        % Generation of Leads and Lagged differences of regressors:
        dx_all = GenLeadLag(DeltaX,lag_comp,lag_comp);
        % Appropriate truncation of all regressors:
        AllUntrunc = [Zs; dx_all];
        AllTrunc = AllUntrunc(:,lag_comp+1:end-lag_comp);
        ytrunc = ys(:,lag_comp+1:end-lag_comp);
        % D-OLS estimation itself:
        DOLSall = (AllTrunc'\ytrunc')';
        % D-OLS residuals:
        u_dols = ytrunc - DOLSall*AllTrunc;
        act_samp = T-1-2*lag_comp;
        SSR = u_dols*u_dols';
        end;
        
        % Computation of information criteria:
        AIC_table(lag_comp+1) = act_samp*log(SSR/act_samp)+2*(m*(2*lag_comp+2)+2);
        BIC_table(lag_comp+1) = act_samp*log(SSR/act_samp)+log(act_samp)*(m*(2*lag_comp+2)+2);
        
        
    end;

    % Computation of optimal lead and lag order:
    %AIC:
    [Asort,Aindex] = sort(AIC_table,'ascend');
    AIC_lag = Aindex(1,1)-1;
    AIC_lead = Aindex(1,1)-1;
    % BIC:
    [Bsort,Bindex] = sort(BIC_table,'ascend');
    BIC_lag = Bindex(1,1)-1;
    BIC_lead = Bindex(1,1)-1;
    
end; % End of symmet loop: