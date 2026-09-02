function [betaGM, betaiGM, V,Vdirect,alphaGM,deltaGM] = GroupMeanFMOLS(y,x,q,type,kern,band,corrrob)
% Group-Mean estimator (mean of N individual specific FM-OLS estimators)
%----------------------------------------------------------------------- 
%INPUTS:    y...            TxN dependent variable
%           x...            TxN integrated regressor
%           q...            highest integer power of regressor (either 2 or 3)
%           type...         either demean (1) or demean and detrend (2) the variables
%           kern/bandw...   kernel and bandwidth for long-run vovariance estimation
%           corrrob...      either 0 (no) or 1 (yes), specifies whether cross-sectional correlation robust VCV estimator should be used
%-----------------------------------------------------------------------
%OUTPUTS:   betaGM...       qx1 group-mean estimate of beta
%           V...            qxq VCV of beta_hat
%           Vdirect...      qxq VCV of beta_hat, direct version
%           alphaGM...      Nx1 estimated individual specific intercepts
%           deltaGM...      Nx1 estimated individual specific time slopes
%------------------------------------------------------------------------
% KR, April 2021
%------------------------------------------------------------------------
%%
% Determine sample length and number of individuals:
[T,N] = size(y);

v = diff(x);
% now, just consider t=2,...,T for the rest of this procedure:
y = y(2:T,:);
x = x(2:T,:);
T = T-1;

% Preallocating matrices:
Eqn(1,N) = struct;

% Transformation of variables:
if q == 2
[ytilde,xtilde,x2tilde,~] = demean_detrend(y,x,type);
elseif q == 3
[ytilde,xtilde,x2tilde,x3tilde] = demean_detrend(y,x,type);
end

% Preallocation:
betaGM = zeros(q,1); 
V = zeros(q,q);
Vdirect = zeros(q,q);

for i = 1:N % Consider the i-th Equation:
    
    % OLS Estimator
    Eqn(i).ytilde = ytilde(:,i);
    Eqn(i).xtilde = xtilde(:,i);
    Eqn(i).x2tilde = x2tilde(:,i);
    Eqn(i).Xtilde = [Eqn(i).xtilde,Eqn(i).x2tilde];
    if q == 3
        Eqn(i).x3tilde = x3tilde(:,i);
        Eqn(i).Xtilde = [Eqn(i).Xtilde,Eqn(i).x3tilde];
    end
    Eqn(i).betaOLS = Eqn(i).Xtilde\Eqn(i).ytilde;

    % OLS Residuals
        % Note that Frisch-Waught Theorem guarantees that this are the
        % correct residuals.
    Eqn(i).u_hat = Eqn(i).ytilde - Eqn(i).Xtilde*Eqn(i).betaOLS;
    Eqn(i).v = v(:,i);
    Eqn(i).mu = mean(Eqn(i).v);
    
    % Long-Run Covariance Matrix Estimation
        % NOTE: for consistency remove potentially non-zero \mu_i from first difference of integrated regressor

        % Calculate bandwidth: 
        %if isequal(band,'And91')
        %    band = And_HAC91([Eqn(i).u_hat, Eqn(i).v - ones(T,1)*mean(Eqn(i).v)],kern);
        %end
        % Kernel estimaton:
        %[Eqn(i).Omega,Eqn(i).Delta,~] = lr_varmod([Eqn(i).u_hat, Eqn(i).v - ones(T,1)*mean(Eqn(i).v)],kern,band,0);
        %javított rész
        if isequal(band,'And91')
            bandw_i = And_HAC91([Eqn(i).u_hat, Eqn(i).v - ones(T,1)*mean(Eqn(i).v)],kern);
        else
            bandw_i = band;
        end
        [Eqn(i).Omega,Eqn(i).Delta,~] = lr_varmod([Eqn(i).u_hat, Eqn(i).v - ones(T,1)*mean(Eqn(i).v)],kern,bandw_i,0);
 
    % Safe the necessary subblocks of Omega and Delta:
    Eqn(i).Omega_uu = Eqn(i).Omega(1,1);
    Eqn(i).Omega_vv = Eqn(i).Omega(2,2);
    Eqn(i).Omega_uv = Eqn(i).Omega(1,2);
    Eqn(i).Omega_vu = Eqn(i).Omega(2,1);
    Eqn(i).Delta_vv = Eqn(i).Delta(2,2);
    Eqn(i).Delta_vu = Eqn(i).Delta(2,1);
    
    Eqn(i).Delta_vu_plus = Eqn(i).Delta_vu - Eqn(i).Delta_vv*(Eqn(i).Omega_vv\Eqn(i).Omega_vu);
    Eqn(i).Omega_udotv = Eqn(i).Omega_uu - Eqn(i).Omega_uv*(Eqn(i).Omega_vv\Eqn(i).Omega_vu);

    % Fully Modified Corrections:
    Eqn(i).ytilde_plus = Eqn(i).ytilde - Eqn(i).v*(Eqn(i).Omega_vv\Eqn(i).Omega_vu);
    Eqn(i).x = x(:,i);
    Eqn(i).C = Eqn(i).Delta_vu_plus*[T,2*sum(Eqn(i).x)]';
    if q == 3
        Eqn(i).C = [Eqn(i).C;Eqn(i).Delta_vu_plus*3*sum(Eqn(i).x.^2)];
    end
    Eqn(i).betaFM = (Eqn(i).Xtilde'*Eqn(i).Xtilde)\(Eqn(i).Xtilde'*Eqn(i).ytilde_plus - Eqn(i).C);
    
    betaGM = betaGM + Eqn(i).betaFM;
    if corrrob == 0
    
    V = V + Eqn(i).Omega_udotv*(eye(q)/(Eqn(i).Xtilde'*Eqn(i).Xtilde));
    
    % direct VCV:
    if type == 1 % individual specific intercepts, only
        if q == 2
            scaling = diag([T^(-3/2),T^(-5/2)]);
            Vdirect = Vdirect + Eqn(i).Omega_udotv*((diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2])/[1/12,1/12;1/12,4/45])*diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2]));
        elseif q == 3
            scaling = diag([T^(-3/2),T^(-5/2),T^(-7/2)]);
            Vdirect = Vdirect + Eqn(i).Omega_udotv*((diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2,(1/Eqn(i).mu)^3])...
                /[1/12,1/12,3/40;1/12,4/45,1/12;3/40,1/12,9/112])...
                *diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2,(1/Eqn(i).mu)^3])); 
        end
    elseif type == 2 % individual specific intercepts and linear time trends
        if q == 2
            scaling = diag([T^(-1),T^(-5/2)]);
            XX = scaling*(Eqn(i).Xtilde'*Eqn(i).Xtilde)*scaling;
            XX(2,2) = (1/180)*Eqn(i).mu^4;
            Vdirect = Vdirect + Eqn(i).Omega_udotv*(eye(q)/XX);
        elseif q == 3
            scaling = diag([T^(-1),T^(-5/2),T^(-7/2)]);
            XX = scaling*(Eqn(i).Xtilde'*Eqn(i).Xtilde)*scaling;
            XX(2:3,2:3) = diag([Eqn(i).mu^2,Eqn(i).mu^3])*[1/180,1/120;1/120,9/700]*diag([Eqn(i).mu^2,Eqn(i).mu^3]);
            Vdirect = Vdirect + Eqn(i).Omega_udotv*(eye(q)/XX);
        end
    end
    
    end

end

betaiGM = zeros(q,N);
for i = 1:N
    betaiGM(:,i) = Eqn(i).betaFM;
end

betaGM = betaGM/N;

if type == 1 % Estimate individual specific intercepts:
    if q == 2
        alphaGM = mean(y,1)' - [mean(x,1)',mean(x.^2,1)']*betaGM;
    elseif q == 3
        alphaGM = mean(y,1)' - [mean(x,1)',mean(x.^2,1)',mean(x.^3,1)']*betaGM;
    end
        % the i-th element equals (ones(T,1)\(y(:,i) - [x(:,i), (x(:,i)).^2, (x(:,i)).^3]*betaGM))';
        % this is used below!
    deltaGM = zeros(N,1);
elseif type == 2
    temp = zeros(N,2);
    if q == 2
        for i=1:N
            temp(i,:) = ([ones(T,1),(1:T)']\(y(:,i) - [x(:,i), (x(:,i)).^2]*betaGM))';
        end
    elseif q == 3
        for i=1:N
            temp(i,:) = ([ones(T,1),(1:T)']\(y(:,i) - [x(:,i), (x(:,i)).^2, (x(:,i)).^3]*betaGM))';
        end
    end
    alphaGM = temp(:,1);
    deltaGM = temp(:,2);
end

if corrrob == 1
%% cross-sectional dependence robust VCV

% Estimate large Omega matrix (2Nx2N):
u_hat = NaN(T,N);
v_hat = NaN(T,N);
for i=1:N
    u_hat(:,i) = Eqn(i).u_hat;
    v_hat(:,i) = Eqn(i).v - ones(T,1)*mean(Eqn(i).v);
end
    
if isequal(band,'And91')
    band = And_HAC91([u_hat, v_hat],kern);
end
[Omega,~,~] = lr_varmod([u_hat, v_hat],kern,band,0); %2Nx2N

Omega_uu = Omega(1:N,1:N);
Omega_uv = Omega(1:N,(N+1):(2*N));
Omega_vu = Omega((N+1):(2*N),1:N);
Omega_vv = Omega((N+1):(2*N),(N+1):(2*N));
    
for i=1:N
    for j=1:N
        
        % standard
        
        Mii = Eqn(i).Xtilde'*Eqn(i).Xtilde;
        Mij = Eqn(i).Xtilde'*Eqn(j).Xtilde;
        Mjj = Eqn(j).Xtilde'*Eqn(j).Xtilde;
        
        Omega_uidotvi_ujdotvj = Omega_uu(i,j) - (Omega_uv(i,i)/Omega_vv(i,i))*Omega_vu(i,j) - ...
                                (Omega_uv(j,j)/Omega_vv(j,j))*Omega_vu(j,i) + ...
                                (Omega_uv(i,i)/Omega_vv(i,i))*Omega_vv(i,j)*(Omega_vv(j,j)\Omega_vu(j,j));
        
        V = V + Omega_uidotvi_ujdotvj*(Mii\(Mij/Mjj));
        
        % direct
        
        if type == 1 % individual specific intercepts, only
            if q == 2
                scaling = diag([T^(-3/2),T^(-5/2)]);
                Vdirect = Vdirect + Omega_uidotvi_ujdotvj*((diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2])/[1/12,1/12;1/12,4/45])*diag([1/Eqn(j).mu,(1/Eqn(j).mu)^2]));
            elseif q == 3
                scaling = diag([T^(-3/2),T^(-5/2),T^(-7/2)]);
                Vdirect = Vdirect + Omega_uidotvi_ujdotvj*((diag([1/Eqn(i).mu,(1/Eqn(i).mu)^2,(1/Eqn(i).mu)^3])...
                    /[1/12,1/12,3/40;1/12,4/45,1/12;3/40,1/12,9/112])...
                    *diag([1/Eqn(j).mu,(1/Eqn(j).mu)^2,(1/Eqn(j).mu)^3])); 
            end
        elseif type == 2 % individual specific intercepts and linear time trends
            if q == 2
                scaling = diag([T^(-1),T^(-5/2)]);
                Xii = scaling*(Eqn(i).Xtilde'*Eqn(i).Xtilde)*scaling;
                Xij = scaling*(Eqn(i).Xtilde'*Eqn(j).Xtilde)*scaling;
                Xjj = scaling*(Eqn(j).Xtilde'*Eqn(j).Xtilde)*scaling;
                Xii(2,2) = Eqn(i).mu^2*(1/180)*Eqn(i).mu^2;
                Xij(2,2) = Eqn(i).mu^2*(1/180)*Eqn(j).mu^2;
                Xjj(2,2) = Eqn(j).mu^2*(1/180)*Eqn(j).mu^2;
                Vdirect = Vdirect + Omega_uidotvi_ujdotvj*(Xii\(Xij/Xjj));
            elseif q == 3
                scaling = diag([T^(-1),T^(-5/2),T^(-7/2)]);
                Xii = scaling*(Eqn(i).Xtilde'*Eqn(i).Xtilde)*scaling;
                Xij = scaling*(Eqn(i).Xtilde'*Eqn(j).Xtilde)*scaling;
                Xjj = scaling*(Eqn(j).Xtilde'*Eqn(j).Xtilde)*scaling;
                Xii(2:3,2:3) = diag([Eqn(i).mu^2,Eqn(i).mu^3])*[1/180,1/120;1/120,9/700]*diag([Eqn(i).mu^2,Eqn(i).mu^3]);
                Xij(2:3,2:3) = diag([Eqn(i).mu^2,Eqn(i).mu^3])*[1/180,1/120;1/120,9/700]*diag([Eqn(j).mu^2,Eqn(j).mu^3]);
                Xjj(2:3,2:3) = diag([Eqn(j).mu^2,Eqn(j).mu^3])*[1/180,1/120;1/120,9/700]*diag([Eqn(j).mu^2,Eqn(j).mu^3]);
                Vdirect = Vdirect + Omega_uidotvi_ujdotvj*(Xii\(Xij/Xjj));
            end
        end  
    end
end

end

% Final VCV scaling:
V = V/(N^2);
Vdirect = scaling*(Vdirect/(N^2))*scaling;


end

