%-------------------------------------------------------------------
%CTPUTests
% function[PO_c,PO_t,C_val_all,t_val_all,coeff_dec,t_dec] .. 
%                  =PhilOul(uhat,d,n,kern,band,deme,alpha)
%
% Function to compute the Phillips and Ouliaris - Phillips and Perron
% type tests (coefficient and t-statistic version).
%
% Input:    uhat ... Residuals from cointegrating OLS regression, R^{T
%                    times 1}
%              d ... Specification deterministics (-1=nothing, 0=intercept,
%                    1=intercept and linear trend)
%              n ... Number of (integrated and not coint.) regressors
%           kern ... Choice of kernel in computations
%                    (tr,ba,pa,bo,da,qs)
%           band ... Choice of bandwidth
%           deme ... Demeaning of residuals in LR-Var. computation
%                    (0=no,1=yes)
%          alpha ... Significance level of test: 0.01, 0.025, 0.05, 0.075, 0.1
%
% Output:       PO_c ... PO Coefficient test statistic. (Z_alpha)
%               PO_t ... PO t-test statistic (Z_t)
%          C_val_all ... Vector (1 times 5) with critical values coeff
%                        test.
%          t_val_all ... Vector (1 times 5) with critical values t-
%                        test.
%          coeff_dec ... Test decision (1=rejection, 0=no rej.) for
%                        Coeff.test at level alpha
%              t_dec ... Test decision (1=rejection, 0=no rej.) for
%                        t-test at level alpha
% External procedures: lr_var, lr_weights
% MW; September 2009
%--------------------------------------------------------------------
function[PO_c,PO_t,C_val_all,t_val_all,coeff_dec,t_dec]... 
                  =PhilOul(uhat,d,n,kern,band,deme,alpha)

    
    [T,shouldbe1] = size(uhat);
    deter = d;

    Yvec = uhat(2:T);
    Xvec = uhat(1:T-1);
    %ahat = inv(Xvec'*Xvec)*Xvec'*Yvec;
    ahat = Xvec\Yvec;
    Uvec = Yvec - Xvec*ahat;
    
    % Variance and long-run variances:
    if isequal(band,'And91');
    bandw = And_HAC91(Uvec,kern);
    elseif isequal(band,'NW');
    bandw = bwNW(Uvec,kern,0,[]);
    else
    bandw = band;
    end
    
  
    [Omega,Delta,Sigma] = lr_var(Uvec,kern,bandw,deme);
    
    %t-value:
    t_ahat = (ahat-1)/sqrt(Sigma*inv(Xvec'*Xvec));
    
    % Phillips-Ouliaris Z_alpha and Z_t tests (always on deter=-1
    % regression)
    PO_c = T*(ahat-1) - 0.5*(Omega-Sigma)/(T^(-2)*Xvec'*Xvec);
    PO_t = sqrt(Sigma)/sqrt(Omega)*t_ahat - 0.5*(Omega-Sigma)*inv(sqrt(Omega)*sqrt(T^(-2)*Xvec'*Xvec));
   
    
    % Tables with critical values (taken from Phillips and Ouliaris,
    % 189--190)
    % Format of Tables: 0.01, 0.025, 0.05, 0.075, 0.10, as columns
    %                   n = 1,2,3,4,5 as rows.
    
    
    if deter == -1,     % No deterministic components:
          tab_coeff=[-22.8291 -18.8833 -15.6377 -13.8123 -12.5438;
                     -29.2688 -25.2101 -21.4833 -19.6142 -18.1785;
                     -36.1619 -31.5432 -27.8526 -25.5236 -23.9225;
                     -42.8724 -37.4769 -33.4784 -30.9288 -27.3952;
                     -48.5240 -42.5473 -38.0934 -35.5142 -32.2654];
                 
              tab_t=[-3.3865 -3.0547 -2.7619 -2.5822 -2.4505;
                     -3.8395 -3.5484 -3.2667 -3.1105 -2.9873;
                     -4.3038 -3.9895 -3.7371 -3.5716 -3.4446;
                     -4.6720 -4.3798 -4.1261 -3.9482 -3.8068;
                     -4.9897 -4.6676 -4.3999 -4.2521 -4.1416];
        
    elseif deter == 0,  % Intercept:
          tab_coeff=[-28.3218 -23.8084 -20.4935 -18.4836 -17.0390;
                     -34.1686 -29.7354 -26.0943 -23.8739 -22.1948;
                     -41.1348 -35.7116 -32.0615 -29.5083 -27.5846;
                     -47.5118 -41.6431 -37.1508 -34.7110 -32.7382;
                     -52.1723 -46.5344 -41.9388 -39.1100 -37.0074];
                 
              tab_t=[-3.9618 -3.6420 -3.3654 -3.1982 -3.0657;
                     -4.3078 -4.0217 -3.7675 -3.5846 -3.4494;
                     -4.7325 -4.3747 -4.1121 -3.9560 -3.8329;
                     -5.0728 -4.7075 -4.4542 -4.2883 -4.1565;
                     -5.2812 -4.9809 -4.7101 -4.5553 -4.4309];
        
    elseif deter == 1,  % Intercept and linear trend:
          tab_coeff=[-35.4185 -30.8451 -27.0866 -24.7530 -23.1915;
                     -40.3427 -36.1121 -32.2231 -29.7331 -27.7803;
                     -47.3590 -42.5998 -37.7304 -34.9951 -33.1637;
                     -53.6142 -47.1068 -42.4593 -39.7286 -37.7368;
                     -58.1615 -52.4874 -47.3830 -44.5074 -42.3231];
                 
              tab_t=[-4.3628 -4.0722 -3.8000 -3.6467 -3.5184;
                     -4.6451 -4.3854 -4.1567 -3.9754 -3.8429;
                     -5.0433 -4.7699 -4.4895 -4.3198 -4.1950;
                     -5.3576 -5.0180 -4.7423 -4.5837 -4.4625;
                     -5.5849 -5.3056 -5.0282 -4.8695 -4.7311];
        
    end;
    
    % Computation of test decisions:
    
    % Selection of row corresponding to the number of regressors for all
    % significance levels:
    C_val_all = tab_coeff(n,:);
    t_val_all = tab_t(n,:);
    
    % Computation of test decision for chosen significance level:
    if alpha == 0.01,
        a_ind = 1;
    elseif alpha == 0.025,
        a_ind = 2;
    elseif alpha == 0.05,
        a_ind = 3;
    elseif alpha == 0.075,
        a_ind = 4;
    elseif alpha == 0.1,
        a_ind = 5;
    end;
    
    C_val_alpha = tab_coeff(n,a_ind);
    t_val_alpha = tab_t(n,a_ind);
    
    if PO_c <= C_val_alpha,
        coeff_dec = 1;
    else
        coeff_dec = 0;
    end;
    
    if PO_t <= t_val_alpha,
        t_dec = 1;
    else
        t_dec = 0;
    end;
    
    