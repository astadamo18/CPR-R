%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CTPUTests
% [CTstat,CTdecis] = CT_test(uplus,omega,d,m,p,alphavec);
% CT_Test: 
%  KPSS-Shin type test for CPRs
%
%  Input: uplus...FM-OLS residuals, T \times 1
%         omega...Estimated Long run variance \omega_u.v
%         d.......Deterministics (-1=none, 0=interc,1=interc+trend)
%         m.......Number of integrated regressors (critical values up 4 regressors)
%         p.......Highest power (critical values up to power 4)
%         alphavec...vector in R^{numAlpha \times 1} percentiles
%
% Output: CTstat....Test statistic
%         CTdecis...Vector in R^{numAlpha \times 1} with rejction
%                   frequencies (0=no rej., 1=rej.)
%
% Note: Percentile vector: [0.01 0.025 0.05 0.1 0.5 0.9 0.95 0.975 0.99]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[CTstat,CTdecis] = CT_test(uplus,omega,d,m,p,alphavec);

% Sample size:
[Tmin1,ddd] = size(uplus);
[numA,dddd] = size(alphavec); % numA is number of percentiles:
 
% Scaling by original length:
T = Tmin1+1;  % Original length.
% Scaling by actual length:
T = Tmin1;   


% Computation of Test Statistic:
partsum = (1/sqrt(T))*cumsum(uplus);
CTstat = (1/(T*omega))*partsum'*partsum;

%sumterm = cumsum(uplus);
%squaredterm = sumterm.^2;
%CT2 = (1/T^2*omega)*sum(squaredterm);

%CTstat
%CT2
 
% Loading of critical values:
if d == 0 || d == 1
    filename_critval = sprintf('CTcritval\\CT_d_%d_m_%d_p_%d.mat',d,m,p);
    load(filename_critval);
    eval(sprintf('crittable = CT_d_%d_m_%d_p_%d;',d,m,p));
    eval(sprintf('clear CT_d_%d_m_%d_p_%d',d,m,p));
elseif d == -1
    filename_critval = sprintf('CTcritval\\CT_d_MIN%d_m%d_p_%d.mat',abs(d),m,p);
    load(filename_critval);
    eval(sprintf('crittable = CT_d_MIN%d_m%d_p_%d;',abs(d),m,p));
    eval(sprintf('clear CT_d_MIN%d_m%d_p_%d',abs(d),m,p));
end


% Test decisions at alpha levels:
for jj = 1:numA;
    
    if alphavec(jj,1) == 0.1;
            critval = crittable(1,6);
        elseif alphavec(jj,1) == 0.05;
            critval = crittable(1,7);
        elseif alphavec(jj,1) == 0.01;
            critval = crittable(1,9);
        end;
    
CTdecis(jj,:) = (CTstat > critval);
        
end;        

