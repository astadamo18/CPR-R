%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CTPUTests
% [PUstat,PUdecis] = PU_test(y,x,d,m,p,kern,band,deme,alphavec);
% PU_Test:
% Null hypothesis is: NO COINTEGRATION 
%
%  Input: y ... dependent variable in R^{T \times 1}
%         x ... Single integrated regressor in R^{T \times m}
%               Regressors 1,...,m-1 with power one
%               Regressor m with powers 'orders'
%         d.......Deterministics (-1=none, 0=interc,1=interc+trend)
%         m.......Number of integrated regressors
%         p.......Highest power
%         orders..If scalar: Highest power of Regressor x_m, then powers 1,...,orders are included
%                 If cell array: powers orders{1} are included for the m-th regressor
%           Note: If cell array, then orders must consist one row-matrix (in order to be in line
%                 with the function GenVarPolyTerms)
%         kern....Kernel function
%         band....Bandwdith choice (Number or NW,And91,AM92)
%         deme....Residual demeaning: choose 0
%         alphavec...vector in R^{numAlpha \times 1} percentiles
%
% Output: PUstat....Test statistic
%         PUdecis...Vector in R^{numAlpha \times 1} with rejction
%                   frequencies (0=no rej., 1=rej.)
%
% Required files: GenVarPolyTerms.m, And_HAC91.m, bwNW.m, AndMon_HAC92.m,
% lr_var.m, trimr.m
%
% Note: Percentile vector: [0.01 0.025 0.05 0.1 0.5 0.9 0.95 0.975 0.99]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[PUstat,PUdecis] = PU_test(y,x,d,m,orders,kern,band,deme,alphavec);

% Sample size:
[T,ddd] = size(y);
[numA,dddd] = size(alphavec); % numA is number of percentiles:


% We consider the case with intercept and linear trend:
const = ones(T,1);
trend = cumsum(const);
d2 = d+2;
switch d2
    case 1 
       deter = [];
    case 2 
       deter = const;
    case 3
        deter = [const trend];
end

% Stacking y and x:
z = [y x];

% Generation of powers of regressor:

% Take the last column of x with powers
xpower = x(:,m);

if isnumeric(orders) && orders == 1;
    powerreg = [];
else reg = GenVarPolyTerms(xpower,orders);
    powerreg = reg.X;
end;

% Computation of residuals u-hat:

if isnumeric(orders) && orders == 1;
    regmat4uhat = [deter x];
else
%powerreg_ = powerreg(:,2:end);
x_ = x(:,1:(m-1));
% Remark: The 1st column of powerreg is already the last column in x, 
% so cut last column of x
regmat4uhat = [deter x_ powerreg];
end 

coeff4uhat = regmat4uhat\y;
uhat = y - regmat4uhat*coeff4uhat;


% Computation of long-run variance using the VAR(1) residuals:
% var of order 1 including the deterministics:

depvar = z(2:end,:)';
lagvar = lag(z);
indepvar = [deter(2:end,:)'; lagvar(2:end,:)'];

%size(indepvar)
%size(depvar)
varcoeff = (indepvar'\depvar')';

varresid = depvar - varcoeff*indepvar;
varresid = varresid';

% (1) Constructing LR variance Estimators
if isequal(band,'And91');
    bandw = And_HAC91(varresid,kern);
elseif isequal(band,'NW');
    bandw = bwNW(varresid,kern,0,[]);
else
    bandw = band;
end

if isequal(band,'AM92');
    [Lr,Dr,Sr] = AndMon_HAC92(varresid,kern,1,1,deme);
else
    [Lr,Dr,Sr] = lr_var(varresid,kern, bandw, deme);
end
%disp('From FM_CPR'),Lr, Dr, Sr
omega_udotv = Lr(1,1)-Lr(1,2:end)*inv(Lr(2:end,2:end))*Lr(2:end,1);

%Lr
%omega_udotv
% Long-run variance computation:
%[Omega_u,Delta_u,Sigma_u] = lr_var(uhat,kern,band,deme);



% Computation of Test Statistic:
PUstat = omega_udotv/(T^(-2)*uhat'*uhat);

% Loading of critical values:
if isnumeric(orders)
if d == 0 || d == 1
    filename_critval = sprintf('PUcritval\\PU_d_%d_m_%d_p_%d.mat',d,m,orders);
    load(filename_critval);
    eval(sprintf('crittable = PU_d_%d_m_%d_p_%d;',d,m,orders));
    eval(sprintf('clear PU_d_%d_m_%d_p_%d',d,m,orders));
elseif d == -1
    filename_critval = sprintf('PUcritval\\PU_d_MIN%d_m%d_p_%d.mat',abs(d),m,orders);
    load(filename_critval);
    eval(sprintf('crittable = PU_d_MIN%d_m%d_p_%d;',abs(d),m,orders));
    eval(sprintf('clear PU_d_MIN%d_m%d_p_%d',abs(d),m,orders));
end


% Test decisions at alpha levels:
for jj = 1:numA;
    
    if alphavec(jj,1) == 0.10;
            critval = crittable(1,6);
        elseif alphavec(jj,1) == 0.05;
            critval = crittable(1,7);
        elseif alphavec(jj,1) == 0.01;
            critval = crittable(1,9);
    end
    
PUdecis(jj,:) = (PUstat > critval);
end
else PUdecis = [];
end %end isnumeric(orders)        
end        

