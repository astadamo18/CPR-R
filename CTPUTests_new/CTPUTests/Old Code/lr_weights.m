%----------------------------------------------------------
% function [w,upper] = lr_weights(T,kern,band)
%
% This function computes the weights corresponding to some
% kernel functions. Only the weights for positive arguments 
% are computed.
%       
% Input:   T     ... Time series dimension of data
%          kern  ... Specifies chosen kernel function:
%                tr ... Truncated
%                ba ... Bartlett
%                pa ... Parzen
%                bo ... Bohman
%                da ... Daniell
%                qs ... Quadratic Spectral
%          band  ... specifies the bandwidth chosen, integer 
%                    in the range 1,...,T                                   
%
% Output:  w     ... Vector of size T-1 containing the weights.
%          upper ... Index to largest non-zero entry in w (to only sum
%                    non-zero terms) in lr_var procedure
% Remark:   USAGE OF BANDWIDTH BAND IS now as IN (I BELIEVE) OTHER CODES AUGMENTED BY 1 TO 
%           GET BAND+1 AS BANDWIDTH BELOW: NOTATIONAL/CONVENTIONAL ISSUE
%
% External function: none
%
% TV and MW, August 2009.
%----------------------------------------------------------
function[w,upper] = lr_weights(T,kern,band)

% Initialization:
w = zeros(T-1,1);
%M = band+1;   %number of terms up to M in truncation kernels.
M = band;
if kern == 'tr',               % Truncated
    for j = 1:M;
        w(j) = 1;
    end;
    upper = min(M,T-1);
elseif kern == 'ba',           % Bartlett
    for j = 1:ceil(M)-1;
        w(j) = 1-j/(M+0);
    end;
    upper = ceil(M)-1;
elseif kern == 'pa',           % Parzen 
    for j = 1:floor(M/2);
        jj = j/M;
        w(j) = 1 - 6*jj^2 + 6*jj^3;
    end;
    for j = floor(M/2)+1:M,
        jj = j/M;
        w(j) = 2*(1-jj)^3;
    end;
    upper = ceil(M) - 1;
elseif kern == 'bo',            % Bohman
    for j =1:ceil(M)-1;
        jj  = j/M;
        w(j) = (1-jj)*cos(pi*jj)+sin(pi*jj)/pi;
    end;
    upper = ceil(M)-1;
elseif kern == 'da',           % Daniell
    for j = 1:T-1;
        w(j) = sin(pi*j/M)/(pi*j/M);
    end;
    upper = T-1;
elseif kern == 'qs',           % Quadratic spectral
    sc = (6*pi)/5;
    for j = 1:T-1;
        jj = j/M;
        w(j) = 25/(12*(pi^2)*(jj^2))*(sin(sc*jj)/(sc*jj)-cos(sc*jj));
    end;
    upper = T-1;
end;

