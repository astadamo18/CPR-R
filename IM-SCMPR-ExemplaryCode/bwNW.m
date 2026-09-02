function result = bwNW(v,kern,inter,weights)
%------------------------------------------------
% Function to compute the Newey West data-depended automatic bandwidth 
% according to Newey West (1994)
%
% In our situation: hat(h)_t = [hat(u)_t ; v_t]
%
% v...     (Txm)-matrix (in our case this will be (û v), where v_t= x_t - x_(t-1) 
% kern...  kernel chosen ('ba','pa','qs')
% inter... Dummy variable: 0... no intercept
%                          1... intercept (normally inter should be set to 0 in our cases)
% weights..(mx1)-vector of weights
%          Default: (1,1,...1)' when there is no intercept,
%                   (0,1,...1)' when there is an intercept
%-----------------------------------------------

T = size(v,1);
m = size(v,2);

% Now transpose v to be in line with first version ( now vmat \in (mxT) )

vmat = v';

% Weightening-vector w (mx1) : 
% Set the vector (1,...1)' as the default when there is no intercept,
% (0,1,...1)' when there is an intercept

if isempty(weights)
    weights = ones(m,1);
    if inter == 1
        weights(1) = 0;
    end
end

% Select the parameter "n" (lag truncation parameter) of NeweyWest (1994) 
% depending on the kernel

switch kern
    case 'ba'
        npower = 2/9;
    case 'pa'
        npower = 4/25;
    case 'qs' 
        npower = 2/25;
end
n = floor(4 * (T/100)^npower);

% Computing Vector containing sigmahat^2_j for j in 0:n (3.9) NW(1994)
vmatw = weights'*vmat;
sigma = zeros(n+1,1);
for j = 1:n+1
    sigma(j) = T^(-1) * sum(vmatw(:,j:T)*vmatw(:,1:(T-j+1))');
end
% Computing shat^(q) for q=0,1,2  (3.10) NW(1994)
s0 = sigma(1) + 2*sum(sigma(2:end));
s1 = 2 * sum((1:n) * sigma(2:end));
s2 = 2 * sum((1:n).^2 * sigma(2:end));

% Finally compute gamma and bandwith according to Lemma 1 (3.7) using 
% the values for q and c_gamma from Table I B (1)&(2) NW (1994)

if kern == 'ba'
    q = 1;
else q = 2;
end
Tpower = 1/(2*q + 1);

switch kern
    case 'ba' 
        gamma = ( 1.1447 * ((s1/s0)^2)^Tpower );
    case 'pa'
        gamma = ( 2.6614 * ((s2/s0)^2)^Tpower );
    case 'qs' 
        gamma = ( 1.3221 * ((s2/s0)^2)^Tpower );
end

bandw = gamma * T^(Tpower);
bandw = min(bandw,n);

% Take integer Bandwitdh 
% bandw = floor(bandw)+1   ;

result = bandw;
return 