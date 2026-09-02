function result = PanelEKC_two_eff(y,x,q,kern,band)
%------------------------------------------------------------------------
% This function is based on de Jong and Wagner (2016). We consider here the
% regressor with its powers up to power 3 and a model where idividual 
% and time effects are present.
%
% INPUTS     y...   left-hand-side (dependent) variable of size TxN
%
%            x...   right-hand-side var. of size TxN (note that we only consider 
%                   the case of k=1 regressor entering with its powers up
%                   to power 3.
%   kern/bandw...   kernel and bandwidth for long-run var estimation
%    
%-------------------------------------------------------------------------
%
% OUTPUTS:  beta_lsdv...  The LSDV-Estimator of beta dim=(qx1) 
%           u_hat...      LSDV residuals (N*T x 1)
%           beta_Mod...   The Modified OLS Estimator (qx1)
%           beta_FM...    The Fully Modified OLS Estimator (qx1)
%           VCV_Mod...    VCV-Matrix for beta_Mod
%           VCV_FM...     VCV-Matrix for beta_FM
%           VCV_FM_std... VCV-Matrix for beta_FM in the standard setting
%                         (i.e. Omega_i = Omega, for all i=1,...,N, 
%                         and Omega non-stochastic)
%           EqnInfo...    A structure array containing all variance
%                         estimators and correction terms (for each i)
%------------------------------------------------------------------------
% KR, March 2017
%------------------------------------------------------------------------
%%
switch q
case 3
% Determine time length and number of individuals
[T,N] = size(y);

% Preallocating correction terms:
Eqn(1,N) = struct;
Sum_Lr = zeros(2,2);
Sum_Dr = zeros(2,2);
Sum_Mi = zeros(q,1);
Sum_FM_cor = zeros(q,1);
Sum_DMD = zeros(q,q);
Sum_ODMD = zeros(q,q);
Sum_ODQD = zeros(q,q);
Sum_Sigma_13 = zeros(q,q);
Sum_OuuOvv = 0;
Sum_Omega_udotv = 0;
Sum_OudotvOvv = 0;

% Demeaning:
[y_check,x_check,xquad_check,xcub_check,~,~,~,~] = simpledemean(y,x,'twoway',q);

% Reshape demeaned TxN matrices into NTx1 vectors
ycheck = vec(y_check);
xcheck = vec(x_check);
x2check = vec(xquad_check);
x3check = vec(xcub_check);

%% OLS (LSDV) Estimation:
Xcheck = [xcheck, x2check ,x3check];% NTxq
XX_check = Xcheck'*Xcheck;
invXX_check = inv(XX_check);
beta_lsdv = Xcheck\ycheck;

u_hat = ycheck-Xcheck*beta_lsdv;

%sigma2 = u_hat'*u_hat/(N*T-(N+T+q));  %????? "fallen nicht parameter weg", da wir demeaned haben?

%% Define needed Matrices:

% First differences of x-matrix:
vt = [x(1,:);diff(x)]; % since x_i0=0

%Weighting matrix:
GT = diag([T^(-1),T^(-1.5),T^(-2)]); %qxq
%Matrix M:
M = [1/6,0,3/8;0,5/12,0;3/8,0,39/20]; %qxq
%Matrix Q:
Q = [1/3,0,9/10;0,59/60,0;9/10,0,101/20]; %qxq

%% Compute individual-wise variables:
for i=1:N
    
Eqn(i).x = x(:,i); 
Eqn(i).u_hat = u_hat(((i-1)*T+1):(i*T)); 
Eqn(i).vt = vt(:,i);
Eqn(i).X_check = Xcheck(((i-1)*T+1):(i*T),:);
   
% Unit-wise bandwidth computation:
if isequal(band,'And91');
    bandw = And_HAC91([Eqn(i).u_hat Eqn(i).vt],kern);
elseif isequal(band,'NW');
    bandw = bwNW([Eqn(i).u_hat Eqn(i).vt],kern,0,[]);
else
    bandw = band;
end

% Long-run, Half-long-run % contemp. variance estimation (indi-wise):
[Eqn(i).Lr,Eqn(i).Dr,Eqn(i).Sr] = lr_varmod([Eqn(i).u_hat, Eqn(i).vt],kern,bandw,0);
    % note: Eqn(i).Lr = Omega_i, Eqn(i).Dr = Delta_i, Eqn(i).Sr = Sigma_i
    
    % Compute sum of Omega_uui*Omega_vvi and Omega_vvi^2, then compute
    % averages outside the loop:
    Sum_OuuOvv = Sum_OuuOvv + Eqn(i).Lr(1,1)*Eqn(i).Lr(2,2); 
    
    % Compute sum of LR_VARs, then compute the averages outside the loop:
    Sum_Lr = Sum_Lr + Eqn(i).Lr;
    Sum_Dr = Sum_Dr + Eqn(i).Dr;    
    
% Correction for Modified OLS:
    % Compute M_i:
    Eqn(i).M = [T;2*sum(Eqn(i).x);3*sum((Eqn(i).x).^2)];
    % Compute the sum of the M_i (stepwise)
    Sum_Mi = Sum_Mi + Eqn(i).M;   
    
% For Inference purposes:
    Eqn(i).D = diag([Eqn(i).Lr(2,2)^(0.5),Eqn(i).Lr(2,2),Eqn(i).Lr(2,2)^(1.5)]);
    Sum_DMD = Sum_DMD + Eqn(i).D*M*Eqn(i).D;
    
    Eqn(i).Omega_udotv = Eqn(i).Lr(1,1)-(Eqn(i).Lr(2,2)^(-1))*Eqn(i).Lr(2,1)^2;
    Sum_Omega_udotv = Sum_Omega_udotv + Eqn(i).Omega_udotv;
    Sum_OudotvOvv = Sum_OudotvOvv + Eqn(i).Omega_udotv*Eqn(i).Lr(2,2);
    
    Sum_ODMD = Sum_ODMD + Eqn(i).Omega_udotv*(Eqn(i).D*M*Eqn(i).D);
    
    % Compute Omega_{uvi}^2Omega_{vvi}^{-1}:
    Eqn(i).Ouv2vv = (Eqn(i).Lr(1,2)^2)*(Eqn(i).Lr(2,2))^(-1);
    Sum_ODQD = Sum_ODQD + Eqn(i).Ouv2vv*(Eqn(i).D*Q*Eqn(i).D);
    
    Sum_Sigma_13 = Sum_Sigma_13 + [0.25*Eqn(i).Lr(1,2)^2,0,0.5*Eqn(i).Lr(2,2)*Eqn(i).Lr(1,2)^2;...
                                   0,0,0;...
                                   0.5*Eqn(i).Lr(2,2)*Eqn(i).Lr(1,2)^2,0,(Eqn(i).Lr(2,2)^2)*(Eqn(i).Lr(1,2)^2)];
 
end % i-loop   
  
%% Compute Correction Terms 

% Compute averages over indivdual-specific LR-VARs:
Lr_mean = Sum_Lr/N;
Dr_mean = Sum_Dr/N;

OuuOvv_mean = Sum_OuuOvv/N;

% Correction for Modified OLS:
    % Compute Sum i=1,...,N \tilde{C}^*_i:
    Sum_Ci_tilde_star = Dr_mean(2,1)*Sum_Mi + N*[-0.5*T*Lr_mean(1,2);0;-(T^2)*Lr_mean(2,2)*Lr_mean(1,2)];

% Correction for FM OLS:
    Dr_vu_plus = Dr_mean(2,1) - Dr_mean(2,2)*((Lr_mean(2,2))^(-1))*Lr_mean(2,1);
    
    for i=1:N
        Eqn(i).C_plus_star = Dr_vu_plus*Eqn(i).M; 
        Eqn(i).y_checkplus_star = y_check(1:end,i) - Lr_mean(1,2)*(Lr_mean(2,2)^(-1))*Eqn(i).vt;
        Eqn(i).FM_cor = Eqn(i).X_check'*Eqn(i).y_checkplus_star - Eqn(i).C_plus_star;
        Sum_FM_cor = Sum_FM_cor + Eqn(i).FM_cor;
    end % i-loop    
    
%% Estimators of beta

% Modified OLS Estimator:
beta_Mod = beta_lsdv - invXX_check*Sum_Ci_tilde_star;
    
% Fully Modified OLS Estimator:
beta_FM = invXX_check*Sum_FM_cor;  
    

%% Estimators of the VCV-Matrices

% VCV OLS (naive):
%VCV_OLSnaive = sigma2*invXX_check;

% VCV Mod OLS:
    V1_hat = Sum_DMD/N;
    V2_hat = V1_hat - diag([0,(Lr_mean(2,2)^2)/12,0]);
    invV2_hat = inv(V2_hat);
    
    % Estimator for Sigma_11, Sigma_12, Sigma_13, Sigma_1 and Sigma2
    Sigma11 = Sum_ODMD/N;
    Sigma12 = Sum_ODQD/N;
    Sigma13 = Sum_Sigma_13/N;
    Sigma1 = Sigma11 + Sigma12 - Sigma13;
                  
Sigma2 = Sigma1 - diag([0,OuuOvv_mean*Lr_mean(2,2)/6,0]) +...
                      diag([0,(Lr_mean(1,1)*Lr_mean(2,2)^2)/12,0]);      

 VCV_Mod = (N^(-1))*GT*(invV2_hat*Sigma2*invV2_hat)*GT;

% VCV FM OLS:
    OudotvOvv_mean = Sum_OudotvOvv/N;
    Omega_udotv_mean = Sum_Omega_udotv/N;
    Sigma2plus = Sigma11 - diag([0,Lr_mean(2,2)*OudotvOvv_mean/6,0]) +... %note: Sigma1plus=Sigma11
                           diag([0,Omega_udotv_mean*(Lr_mean(2,2)^2)/12,0]); 
               
VCV_FM = (N^(-1))*GT*(invV2_hat*Sigma2plus*invV2_hat)*GT;

% standard FM-OLS Version
    % (i.e. Omega_i = Omega, for all i=1,...,N, and Omega non-stochastic)
VCV_FM_std = Omega_udotv_mean * invXX_check;

%% Save Results
% result.V2_hat = V2_hat;
% result.Sigma2 = Sigma2;
% result.Sigma2plus = Sigma2plus;
% result.u_hat = u_hat;

result.ycheck = ycheck;
result.xcheck = xcheck;
result.beta_lsdv = beta_lsdv;
result.beta_Mod = beta_Mod;
result.beta_FM = beta_FM;
result.EqnInfo = Eqn;
%result.VCV_OLSnaive = VCV_OLSnaive;
result.V2_hat = V2_hat;
result.Sigma2 = Sigma2;
result.VCV_Mod = VCV_Mod;
result.VCV_FM = VCV_FM;
result.VCV_FM_std = VCV_FM_std;
    

case 2
% Determine time length and number of individuals
[T,N] = size(y);

% Preallocating correction terms:
Eqn(1,N) = struct;
Sum_Lr = zeros(2,2);
Sum_Dr = zeros(2,2);
Sum_Mi = zeros(q,1);
Sum_FM_cor = zeros(q,1);
Sum_DMD = zeros(q,q);
Sum_ODMD = zeros(q,q);
Sum_ODQD = zeros(q,q);
Sum_Sigma_13 = zeros(q,q);
Sum_OuuOvv = 0;
Sum_Omega_udotv = 0;
Sum_OudotvOvv = 0;

% Demeaning:
[y_check,x_check,xquad_check,~,~,~,~,~] = simpledemean(y,x,'twoway',q);

% Reshape demeaned TxN matrices into NTx1 vectors
ycheck = vec(y_check);
xcheck = vec(x_check);
x2check = vec(xquad_check);

%% OLS (LSDV) Estimation:
Xcheck = [xcheck, x2check];% NTxq
XX_check = Xcheck'*Xcheck;
invXX_check = inv(XX_check);
beta_lsdv = Xcheck\ycheck;

u_hat = ycheck-Xcheck*beta_lsdv;

%sigma2 = u_hat'*u_hat/(N*T-(N+T+q));  %????? "fallen nicht parameter weg", da wir demeaned haben?

%% Define needed Matrices:

% First differences of x-matrix:
vt = [x(1,:);diff(x)]; % since x_i0=0

%Weighting matrix:
GT = diag([T^(-1),T^(-1.5)]); %qxq
%Matrix M:
M = [1/6,0;0,5/12]; %qxq
%Matrix Q:
Q = [1/3,0;0,59/60]; %qxq

%% Compute individual-wise variables:
for i=1:N
    
Eqn(i).x = x(:,i); 
Eqn(i).u_hat = u_hat(((i-1)*T+1):(i*T)); 
Eqn(i).vt = vt(:,i);
Eqn(i).X_check = Xcheck(((i-1)*T+1):(i*T),:);
   
% Unit-wise bandwidth computation:
if isequal(band,'And91');
    bandw = And_HAC91([Eqn(i).u_hat Eqn(i).vt],kern);
elseif isequal(band,'NW');
    bandw = bwNW([Eqn(i).u_hat Eqn(i).vt],kern,0,[]);
else
    bandw = band;
end

% Long-run, Half-long-run % contemp. variance estimation (indi-wise):
[Eqn(i).Lr,Eqn(i).Dr,Eqn(i).Sr] = lr_varmod([Eqn(i).u_hat, Eqn(i).vt],kern,bandw,0);
    % note: Eqn(i).Lr = Omega_i, Eqn(i).Dr = Delta_i, Eqn(i).Sr = Sigma_i
    
    % Compute sum of Omega_uui*Omega_vvi and Omega_vvi^2, then compute
    % averages outside the loop:
    Sum_OuuOvv = Sum_OuuOvv + Eqn(i).Lr(1,1)*Eqn(i).Lr(2,2); 
    
    % Compute sum of LR_VARs, then compute the averages outside the loop:
    Sum_Lr = Sum_Lr + Eqn(i).Lr;
    Sum_Dr = Sum_Dr + Eqn(i).Dr;    
    
% Correction for Modified OLS:
    % Compute M_i:
    Eqn(i).M = [T;2*sum(Eqn(i).x)];
    % Compute the sum of the M_i (stepwise)
    Sum_Mi = Sum_Mi + Eqn(i).M;   
    
% For Inference purposes:
    Eqn(i).D = diag([Eqn(i).Lr(2,2)^(0.5),Eqn(i).Lr(2,2)]);
    Sum_DMD = Sum_DMD + Eqn(i).D*M*Eqn(i).D;
    
    Eqn(i).Omega_udotv = Eqn(i).Lr(1,1)-(Eqn(i).Lr(2,2)^(-1))*Eqn(i).Lr(2,1)^2;
    Sum_Omega_udotv = Sum_Omega_udotv + Eqn(i).Omega_udotv;
    Sum_OudotvOvv = Sum_OudotvOvv + Eqn(i).Omega_udotv*Eqn(i).Lr(2,2);
    
    Sum_ODMD = Sum_ODMD + Eqn(i).Omega_udotv*(Eqn(i).D*M*Eqn(i).D);
    
    % Compute Omega_{uvi}^2Omega_{vvi}^{-1}:
    Eqn(i).Ouv2vv = (Eqn(i).Lr(1,2)^2)*(Eqn(i).Lr(2,2))^(-1);
    Sum_ODQD = Sum_ODQD + Eqn(i).Ouv2vv*(Eqn(i).D*Q*Eqn(i).D);
    
    Sum_Sigma_13 = Sum_Sigma_13 + [0.25*Eqn(i).Lr(1,2)^2,0;0,0];
 
end % i-loop   
  
%% Compute Correction Terms 

% Compute averages over indivdual-specific LR-VARs:
Lr_mean = Sum_Lr/N;
Dr_mean = Sum_Dr/N;

OuuOvv_mean = Sum_OuuOvv/N;

% Correction for Modified OLS:
    % Compute Sum i=1,...,N \tilde{C}^*_i:
    Sum_Ci_tilde_star = Dr_mean(2,1)*Sum_Mi + N*[-0.5*T*Lr_mean(1,2);0];

% Correction for FM OLS:
    Dr_vu_plus = Dr_mean(2,1) - Dr_mean(2,2)*((Lr_mean(2,2))^(-1))*Lr_mean(2,1);
    
    for i=1:N
        Eqn(i).C_plus_star = Dr_vu_plus*Eqn(i).M; 
        Eqn(i).y_checkplus_star = y_check(1:end,i) - Lr_mean(1,2)*(Lr_mean(2,2)^(-1))*Eqn(i).vt;
        Eqn(i).FM_cor = Eqn(i).X_check'*Eqn(i).y_checkplus_star - Eqn(i).C_plus_star;
        Sum_FM_cor = Sum_FM_cor + Eqn(i).FM_cor;
    end % i-loop    
    
%% Estimators of beta

% Modified OLS Estimator:
beta_Mod = beta_lsdv - invXX_check*Sum_Ci_tilde_star;
    
% Fully Modified OLS Estimator:
beta_FM = invXX_check*Sum_FM_cor;  
    

%% Estimators of the VCV-Matrices

% VCV OLS (naive):
%VCV_OLSnaive = sigma2*invXX_check;

% VCV Mod OLS:
    V1_hat = Sum_DMD/N;
    V2_hat = V1_hat - diag([0,(Lr_mean(2,2)^2)/12]);
    invV2_hat = inv(V2_hat);
    
    % Estimator for Sigma_11, Sigma_12, Sigma_13, Sigma_1 and Sigma2
    Sigma11 = Sum_ODMD/N;
    Sigma12 = Sum_ODQD/N;
    Sigma13 = Sum_Sigma_13/N;
    Sigma1 = Sigma11 + Sigma12 - Sigma13;
    Sigma2 = Sigma1 - diag([0,OuuOvv_mean*Lr_mean(2,2)/6]) +...
                      diag([0,(Lr_mean(1,1)*Lr_mean(2,2)^2)/12]);
                  
VCV_Mod = (N^(-1))*GT*(invV2_hat*Sigma2*invV2_hat)*GT;

% VCV FM OLS:
    OudotvOvv_mean = Sum_OudotvOvv/N;
    Omega_udotv_mean = Sum_Omega_udotv/N;
    Sigma2plus = Sigma11 - diag([0,Lr_mean(2,2)*OudotvOvv_mean/6]) +... %note: Sigma1plus=Sigma11
                           diag([0,Omega_udotv_mean*(Lr_mean(2,2)^2)/12]); 
               
VCV_FM = (N^(-1))*GT*(invV2_hat*Sigma2plus*invV2_hat)*GT;

% standard FM-OLS Version
    % (i.e. Omega_i = Omega, for all i=1,...,N, and Omega non-stochastic)
VCV_FM_std = Omega_udotv_mean * invXX_check;

%% Save Results
% result.V2_hat = V2_hat;
% result.Sigma2 = Sigma2;
% result.Sigma2plus = Sigma2plus;
% result.u_hat = u_hat;

result.ycheck = ycheck;
result.xcheck = xcheck;
result.beta_lsdv = beta_lsdv;
result.beta_Mod = beta_Mod;
result.beta_FM = beta_FM;
result.EqnInfo = Eqn;
%result.VCV_OLSnaive = VCV_OLSnaive;
result.VCV_Mod = VCV_Mod;
result.VCV_FM = VCV_FM;
result.VCV_FM_std = VCV_FM_std;    
    
    
    
end