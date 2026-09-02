function[ytilde,xtilde,xquad_tilde,xcub_tilde,ymeans,xmeans,x2means,x3means] = simpledemean(y,x,way,q)
% ---------------------------------------------------------------------------
% INPUTS: 
% x : TxN matrix
% y : TxN matrix
% way : Either 'oneway' or 'twoway' demean procedure
% q : highest power of regressor. either 2 or 3.

% OUTPUTS:
% ytilde (TxN) matrix
% xtilde (TxN) matrix
% xquad_tilde (TxN) matrix
% xcub_tilde (TxN) matrix
% ymeans (Nx1) matrix
% xmeans (Nx1) matrix
% x2means (Nx1) matrix
% x3means (Nx1) matrix

    [T,N] = size(x);

if strcmp(way,'oneway')
%mean(x,1) average of each column (1xN)
%ones(T,1)=(1,...,1)'
%ones(T,1)*xmeans is TxN, with the average of the i-th individual in the i-th column
    x2 = x.^2;
    xmeans = mean(x,1); % averages over T for each individual 1:N, (1xN) vector 
    x2means = mean(x2,1);
    ymeans = mean(y,1); 
    
    xQmean = ones(T,1) * xmeans;
    x2Qmean = ones(T,1) * x2means;
    yQmean = ones(T,1) * ymeans;
    
    xtilde = x - xQmean;
    xquad_tilde = x2 - x2Qmean;
    ytilde = y - yQmean;
    
    xmeans = xmeans';
    x2means = x2means';
    ymeans = ymeans';
  
if q == 3
    x3 = x.^3;
    x3means = mean(x3,1);
    x3Qmean = ones(T,1) * x3means;
    xcub_tilde = x3 - x3Qmean;
    x3means = x3means';
else
    xcub_tilde = NaN;
    x3means = NaN;
end

elseif strcmp(way, 'twoway')
    x2 = x.^2;
    xmeans = NaN;
    ymeans = NaN;
    x2means = NaN;
    x3means = NaN;
    
    xtilde = x - ones(T,1) * mean(x,1) - mean(x,2)* ones(N,1)' + sum(sum(x,1),2)/(N*T);
    xquad_tilde = x2 - ones(T,1) * mean(x2,1) - mean(x2,2)* ones(N,1)' + sum(sum(x2,1),2)/(N*T);
    ytilde = y - ones(T,1) * mean(y,1) - mean(y,2)* ones(N,1)' + sum(sum(y,1),2)/(N*T);
    if q == 3
       x3 = x.^3; 
       xcub_tilde = x3 - ones(T,1) * mean(x3,1) - mean(x3,2)* ones(N,1)' + sum(sum(x3,1),2)/(N*T);
    else
        xcub_tilde = NaN;
    end 
else error('Specify one-way or two-way demeaning!')
end
end


  % the following duration is based on the older version of simpledemean!!!

% Simulation_PanelFM_VCVcheck  813.730 seconds. calls simpledemean
% 48000, takes 159.999 seconds

% Simulation_PanelFM_VCVcheck  2548.844 seconds, calls paneldemean 
% 48000 times, takes 2048.773 s

