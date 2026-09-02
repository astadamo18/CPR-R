%------------------------------------------------------------------------------
% function result =  monitoringCPR(x,y,m,method,deter,kern,band,D_opt,sign_lev,detectortype,degree,n)

% Procedure for Monitoring Cointegrating Polynomial Regressions (CPRs).
%
% Input:       y  ... Dependent variable in R^{T \times 1}.
%              x  ... Integrated regressors in R^{T \times k}, without powers.
%          method ... Specifies the method used for modified OLS calculations 
%                     (for cointegrating polynomial regressions):
%                       FM-OLS ('FM_CPR'), IM-OLS ('IM_CPR') or D-OLS ('D_CPR').
%              m  ... Length of calibration period as fraction of the length
%                       of x (between 0.1 and 0.9) or as number of
%                       observations. Default is m = 0.5.
%          deter  ... Specifies determinstics
%                       0 ... only constant
%                       p ... constant + deterministics with highest order
%                       p, if p > 1. Default is deter = 0.
%           kern  ... Specifies chosen kernel function:
%                       tr ... Truncated
%                       ba ... Bartlett
%                       pa ... Parzen
%                       bo ... Bohman
%                       da ... Daniell
%                       qs ... Quadratic Spectral
%                       Default is the Bartlett kernel. Not every kernel
%                       function works with every bandwidth selection rule!
%           band  ... Specifies the bandwidth chosen,
%                       Either a number or Strings('And91','NW').
%                       Default is band='NW'.
%          D_opt  ... Options for the D-CPR calculations.
%                       .n_lead
%                       .n_lag
%                       .kmax ('k4','k12')
%                       .info_crit ('AIC','BIC')
%                       If not used, set D_opt = [].
%       sign_lev  ... Level of significance (between 0.01 and 0.1).
%                       Default is sign_lev=0.05.
%   detectortype  ... Type of detector to be used for the monitoring
%                       procedure (1, 2, 3, 4 or 5). Default is detectortype = 1.
%                       See Knorre, Wagner and Grupe for detector
%                       definitions.
%         degree  ... Degree of the last regressor (only the last regressor
%                       has higher order powers). Default is 1.
%              n  ... Window length for detectors 4 and 5. Possible values:
%                       n = 0.1, 0.2, 0.3.
%
% Output:
%       result
%         .detector     ... Identifier of the detector  
%         .monitoring_stat  
%                       ... Value of the test statistic at break (if a break 
%                               is detected) or maximum monitoring statistic. 
%         .detection_time
%                       ... Detected time of structural break.
%         .critical_value
%                       ... Critical value of the procedure.
%         .sign_level   ... Significance level used for the procedure.
%         .m_frac       ... Calibration period (fraction).
%         .m_index      ... Calibration period (number of observations).
%         .monitoring_stats_all
%                       ... Values of monitoring statistic for each time
%                               point. In R^{T \times 1}. NaNs for calibration
%                               period.
%         .monitoring_stats_max
%                       ... Maximum of monitoring statistics.
%         .residuals    ... Residuals of modified OLS estimation for 
%                               calculating the monitoring statistic in R^{T-1 \times 1}.
%         .method       ... Modified OLS estimation method ('FM_CPR', 'D_CPR', or 'IM_CPR')
%         .kernel       ... Kernel function (name).
%         .bandwidth_method
%                       ... Bandwidth function (name).
%         .bandwidth    ... Bandwidth parameter.
%         .Omega_udotv  ... Cond. LR-var. based on OLS residuals.
%         .S            ... Residual partial sums of modified OLS estimation  
%                           (calibration period) for calculating the 
%                           monitoring statistic in R^{T-1 \times 1}.
%         .n_lag,n_lead ... Lag and lead used in case of D-CPR estimation.
%
% Details:
% The calibration period can be specified by setting the argument m
% to the number of its last observation.
% The corresponding fraction of the data's length will be calculated
% automatically. Alternatively you can set m directly to the fitting
% fraction value. Attention: The calibration period may become smaller than
% intended: The last observation is calculated as floor(m * T)
% (with T = length of x).
%----------------------------------------------------------
% (Required Files)
%   And_HAC91.m, bwNW.m, DOLS_CPR.m, FM_CPR.m, GenCPRCorrVec.m, GenLeadLag.m, GenPowerReg.m, GenVarPolyTerms.m, 
%   IMOLS_NL.m, lag.m, LagLeadOrdersDOLS_CPR.m, lr_var.m, lr_weights.m, poly_regressors.m, 
%   trimr.m, var.m, Vhat_NEW.m, Z_regressors.m.


% Version 24.03.2020
% MG and FK.
function result = monitoringCPR(x,y,m,method,deter,kern,band,D_opt,sign_lev,detectortype,degree,n)

%% Basics

% Set defaults:
if isempty(m)
    m = 0.5;
end

if isempty(n)
    n = [0.1 0.2 0.3];
end

if isempty(deter)
    deter = 0;
end

if isempty(kern)
    kern = 'ba';
end

if isempty(band)
    band = 'NW';
end

if isempty(sign_lev)
    sign_lev = 0.05;
end

if isempty(detectortype)
    detectortype = 1;
end

if isempty(degree)
    degree = 1;
end
%%
n_vec = [0.1 0.2 0.3]; % Values of n for which critical values have been computed. This vector has to be the same as used for CVs.
                               % Difference to n: n can be a single value and should be an element or subset of n_vec. n_vec represents the settings used for critical values.
    
% Number of integrated regressors and vector that specifies the orders.
% Only the last regressor is allowed to have higher orders.
nreg = size(x, 2);
order = ones(1, nreg);
order(:, nreg) = degree;


% Critival values only for weighting function "1"
weighting_func = 1;
 
% Load tables with critical values:
if strcmp(method,'FM_CPR') || strcmp(method,'D_CPR')
    load(sprintf('CritVals/FM/MonCoint_CV_FM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype));
    crit_vals = eval(sprintf('MonCoint_CV_FM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype));
    eval(sprintf('clear MonCoint_CV_FM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype))
        
elseif strcmp(method,'IM_CPR')
    load(sprintf('CritVals/IM/MonCoint_CV_IM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype));
    crit_vals = eval(sprintf('MonCoint_CV_IM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype));
    eval(sprintf('clear MonCoint_CV_IM_deter_%d_nreg_%d_degree_%d_w_func_%d_detector_%d',deter,nreg,degree,weighting_func,detectortype))
end


T = size(x,1);

% Determine calibration fraction and calibration sample
if m >= 1 && m <= (0.9 * T) && m >= (0.1 * T)
    m_index = floor(m);
elseif m < 1 && m <= 0.9 && m >= 0.1
    mT_rounded = round(m * T * 100000) / 100000; % m*T rounded to three decimal places. Otherwise e.g. for m=0.57 and T=200, floor(m*T)=113 and not 114 as it should. Rounding is recommended by the Mathworks support. 
    m_index = floor(mT_rounded);
else
    error('Error: m out of possible range.');
end

m_frac = m_index / T;

%if m ~= m_frac && m ~= m_index
%    disp(['Warning: Specified parameter m: Observation no. ' num2str(m_index) ...
%        ' is the end of the calibration period'])
%end

%% Deterministics and Regressors:

% Deterministics and exponent for weighting function 
if deter == 0
    determ = ones(T, 1, 1);
    ex = 3;
elseif deter == 1
    determ = [ones(T, 1, 1) bsxfun(@power,cumsum(ones(T, 1, 1), 1),(1:deter))];
    ex = 5;
end

x_k = x(:, nreg);
y_m = y(1:m_index,:);
x_m = x(1:m_index,:);
x_m_k = x_m(:, nreg);
determ_m = determ(1:m_index,:);

% Regressor matrix
if degree == 1
    Z = [determ x];
    Z_m = [determ_m x_m];
else
    Z = [determ x repmat(x_k, [1, degree - 1]).^(repmat(2:degree, [T, 1]))];
    Z_m = [determ_m x_m repmat(x_m_k, [1, degree - 1]).^(repmat(2:degree, [m_index, 1]))];
end

v = diff(x,1,1);
v_m = diff(x_m,1,1);

%% OLS Regression
b_ols = Z_m\y_m;
u_ols = y_m - Z_m*b_ols;

% Bandwidth for long run variance estimation
switch band
    case 'NW'
        bw = bwNW([u_ols(2:end,1) v_m],kern,0,[]);
        bandwidth = 'Newey-West';
    case 'And91'
        bw = And_HAC91([u_ols(2:end,1) v_m],kern);
        bandwidth = 'Andrews';
    otherwise
        bw = band;
        bandwidth = 'Set by user';
end

switch kern
    case 'tr'
        kern_name = 'Truncated';
    case 'ba'
        kern_name = 'Bartlett';
    case 'pa'
        kern_name = 'Parzan';
    case 'da'
        kern_name = 'Daniell';
    case 'qs'
        kern_name = 'Quadratic Spectral';
end
         
% Long run variance estimation
[Lr, ~, ~] = lr_var([u_ols(2:end,1) v_m],kern, bw, 0);
Omega_udotv = Lr(1,1)-Lr(1,2:end)*inv(Lr(2:end,2:end))*Lr(2:end,1);


% FM-CPR
if strcmp(method,'FM_CPR')
    FM_CPR_m = FM_CPR(y_m,x_m,order,[],determ_m,kern,band,0,[]);
    y_plus = y(2:end,:) - v * inv(Lr(2:end,2:end))*Lr(2:end,1);
    X = [x(2:end, :) repmat(x_k(2:end, :), [1, degree - 1]).^(repmat(2:degree, [T-1, 1]))];
    u_plus = y_plus - X * FM_CPR_m.beta_fm - determ(2:end,:) * FM_CPR_m.delta_fm;
    u = [0;u_plus];
    u_out = [NaN;u_plus];
end

% D-CPR
if  strcmp(method,'D_CPR')
    % Lead and lag choices
    k4 = floor(4 * (m_index / 100)^(1/4));
    k12 = floor(12 * (m_index / 100)^(1/4));
    if isempty(D_opt)
        kmax = k4;
        [AIC_lag,AIC_lead,~,~] = LagLeadOrdersDOLS_CPR(y_m',x_m',degree,determ_m',kmax,kmax,0);
        n_lead = AIC_lead;
        n_lag = AIC_lag;
    else
        % Check if both, number of lags (n_lag) and number of leads
        % (n_lead) are provided. If not or empty, determine both.
        if  ~isfield(D_opt,'n_lag')  || ~isfield(D_opt,'n_lead') || isempty(D_opt.n_lag) || isempty(D_opt.n_lead) 
            % Check if maximum number of leads and lags (kmax) is provided
            if ~isfield(D_opt,'kmax')  || isempty(D_opt.kmax)
                kmax = k4;
            elseif strcmp(D_opt.kmax,'k4')
                kmax = k4;
            elseif strcmp(D_opt.kmax,'k12')
                kmax = k12;
            end
            % Check if information criterion is provided
            if ~isfield(D_opt,'info_crit') || isempty(D_opt.info_crit)
                [AIC_lag,AIC_lead,~,~] = LagLeadOrdersDOLS_CPR(y_m',x_m',degree,determ_m',kmax,kmax,0);
                n_lead = AIC_lead;
                n_lag = AIC_lag;
            elseif strcmp(D_opt.info_crit,'AIC')
                [AIC_lag,AIC_lead,~,~] = LagLeadOrdersDOLS_CPR(y_m',x_m',degree,determ_m',kmax,kmax,0);
                n_lead = AIC_lead;
                n_lag = AIC_lag;
            elseif strcmp(D_opt.info_crit,'BIC')
                [~,~,BIC_lag,BIC_lead] = LagLeadOrdersDOLS_CPR(y_m',x_m',degree,determ_m',kmax,kmax,0);
                n_lead = BIC_lead;
                n_lag = BIC_lag;
            end
        % If both lead and leg are provided
        else
            n_lead = D_opt.n_lead;
            n_lag = D_opt.n_lag;
        end
    end
    
    % D-CPR estimation
    result_dols = DOLS_CPR(y_m,x_m,order,determ_m,n_lead,n_lag,kern,bw,0);
    theta_dols = [result_dols.delta_dols' result_dols.beta_dols' result_dols.chi_dols'];

    if n_lag + n_lead == 0
        y_trunc = y;
        u_d = y_trunc - Z * theta_dols';
        u = u_d; 
        u_out = u_d;
    else
        x_delta = diff(x,1,1);
        dx_all = GenLeadLag(x_delta',n_lag,n_lead);
        Zs = Z(2:end,:);
        all_untrunc = [Zs dx_all'];
        T1 = size(all_untrunc,1);
        all_trunc = all_untrunc((n_lag + 1):(T1 - n_lead),:);
        ys = y(2:end,:);
        T2 = size(ys,1);
        y_trunc = ys((n_lag + 1):(T2 - n_lead),:);
        u_d = y_trunc - all_trunc * theta_dols';
        u = [0;zeros(n_lag,1);u_d;zeros(n_lead,1)];
        u_out = [NaN; NaN(n_lag,1);u_d;NaN(n_lead,1)];
    end
end

% IM-CPR. Definition of S.
if strcmp(method,'IM_CPR')
    x_m_k_NL = [];
    if degree >= 2
        x_m_k_NL = bsxfun(@power,x_m_k,(2:degree));
    end
    IM_m = IMOLS_NL(y_m,determ_m, x_m, x_m_k_NL,1,1,1);
    Z_cs = cumsum(Z,1);
    multip = [Z_cs x];
    S = cumsum(y) - multip * IM_m.theta_hat1;
    u_im = diff(S,1,1);
    u_out = [NaN;u_im];
    cv_val = 'IM';
else
    S = cumsum(u);
    cv_val = 'FMD';
end

%% Calculate the monitoring statistics:
S_dev = (S / sqrt(T)).^2 / T;
cumsum_ms = cumsum(S_dev((m_index + 1):T));
sum_1m = sum(S_dev(1:m_index));

s = ((m_index + 1):T) ./ T;

% Weighting function 
w = s'.^ex; 
    
% Vector of window sizes as row vector
n = n(:)';  

% The detectors
if detectortype == 1
    H_s_m = cumsum_ms / Omega_udotv;
elseif detectortype == 2
    H_s_m = (cumsum_ms - sum_1m) / Omega_udotv;
elseif detectortype == 3
    H_s_m = cumsum_ms / sum_1m;
elseif detectortype == 4
    Hsm = zeros(1,size(n,2));
    stats_all = zeros(T,size(n,2));
    for j = 1:size(n,2)
        n_index = floor(n(j) * T);
        cs_S_dev = cumsum(S_dev);
        
        % Calculate number of zeros neccessary  at the beginning of the second cumulative sum
        if m_index-n_index+1 > 0    % We subtract cumsum up to "m_index_i-n_index+1", if the related index is zero or smaller, we have to subtract 0.
            add_zeros = 0;
        else
            add_zeros = (m_index-n_index+1) * (-1) + 1;
        end
        
        H_s_m = (cs_S_dev(m_index+1:end) - [zeros(add_zeros,1); cs_S_dev(max(1,m_index-n_index+1):(end-n_index))]) / Omega_udotv;
        
        Hsm(1,j) = max(abs(H_s_m ./ w));
        stats_all(:,j) = [NaN(m_index, 1); abs(H_s_m ./ w)];
    end
elseif detectortype == 5
    Hsm = zeros(1,size(n,2));
    stats_all = zeros(T,size(n,2));
    for j = 1:size(n,2)
        n_index = floor(n(j) * T);
        cs_S_dev = cumsum(S_dev);
        
        % Calculate number of zeros neccessary  at the beginning of the second cumulative sum
        if m_index-n_index+1 > 0    % We subtract cumsum up to "m_index_i-n_index+1", if the related index is zero or smaller, we have to subtract 0.
            add_zeros = 0;
        else
            add_zeros = (m_index-n_index+1) * (-1) + 1;
        end
        
        cs_S_dev_mov = (cs_S_dev(m_index+1:end) - [zeros(add_zeros,1); cs_S_dev(max(1,m_index-n_index+1):(end-n_index))]); % Without Omega_udotv
        H_s_m = cs_S_dev_mov / sum_1m;
        
        Hsm(1,j) = max(abs(H_s_m ./ w));
        stats_all(:,j) = [NaN(m_index, 1); abs(H_s_m ./ w)];
    end
        
end

if detectortype == 1 || detectortype == 2 || detectortype == 3 
    Hsm = max(abs(H_s_m ./ w));
    stats_all = [NaN(m_index, 1); abs(H_s_m ./ w)];
end

m_vec = round((0.1:0.01:0.9)*100)/100';
cv_which_m = (m_vec == m_frac);

%% Extract specific critical values, rejection decision and detection time

if detectortype == 4 || detectortype == 5
    
    critical_value = zeros(1,size(n,2));
    rej_time = zeros(1,size(n,2));
    Hsm_break_or_max = zeros(1,size(n,2));

    for j = 1:size(n,2)
        n_j = n(j);

        % Selector critical values table dependent on window size n
        cv_which_n = (n_vec == n_j);

        % Check if critical value exists for chosen m and interpolate if not
        if any(cv_which_m)
            crit_val = crit_vals{cv_which_n, 1};
            cv_tab = crit_val(cv_which_m, :);
        else
            crit_val = crit_vals{cv_which_n, 1};
            cv_tab = zeros(1, 4);
            for i = 1:4
                cv_tab(:, i) = spline(m_vec, crit_val(:, i), m_frac);
            end
        end

        % Get or interpolate critical value
        p_table = [0.1 0.05 0.025 0.01];

        if ismember(sign_lev,p_table)
            critical_value(1,j) = cv_tab(ismember(p_table,sign_lev));
        else
            critical_value(1,j) = spline(p_table,cv_tab,sign_lev);
        end

        % Store rejection time or set to infinity
        if any(Hsm(1,j) > critical_value(1,j))
            rej_time(1,j) = find(stats_all(:,j) > critical_value(1,j),1);
            Hsm_break_or_max(1,j) = stats_all(find(stats_all(:,j) > critical_value(1,j),1),j);
        else
            rej_time(1,j) = Inf;
            Hsm_break_or_max(1,j) = Hsm(1,j);
        end

    end % end window loop
    
else % Detectortype = 1,2 or 3 -------------------------------------------%
        % Check if critical value exists for chosen m and interpolate if not
        if any(cv_which_m) 
            cv_tab = crit_vals(cv_which_m, :);
        else
            cv_tab = zeros(1, 4);
            for i = 1:4
                cv_tab(:, i) = spline(m_vec, crit_vals(:, i), m_frac);
            end
        end

        % Get or interpolate critical value
        p_table = [0.1 0.05 0.025 0.01];

        if ismember(sign_lev,p_table)
            critical_value = cv_tab(ismember(p_table,sign_lev));
        else
            critical_value = spline(p_table,cv_tab,sign_lev);
        end

        % Store rejection time or set to infinity
        if any(Hsm > critical_value)
            rej_time = find(stats_all > critical_value,1);
            Hsm_break_or_max = stats_all(find(stats_all > critical_value,1));
        else
            rej_time = Inf;
            Hsm_break_or_max = Hsm;
        end
    
end

    
%% Outputs:
result.detector = detectortype;
if detectortype == 4 || detectortype == 5
   result.window_size = n; 
end
result.monitoring_stat = Hsm_break_or_max;
result.detection_time = rej_time;
result.critical_value = critical_value;
result.sign_level = sign_lev;
result.m_frac = m_frac;
result.m_index = m_index;
result.monitoring_stats_all = stats_all;
result.monitoring_stats_max = Hsm;
result.residuals = u_out;
result.method = method;
result.kernel = kern_name;
result.bandwidth_method = bandwidth;
result.bandwidth = bw;
result.Omega_udotv = Omega_udotv;
result.S = S;
% result.u_ols = u_ols;

if strcmp(method,'D_CPR') 
    result.n_lag = n_lag;
    result.n_lead = n_lead;
end

return;
