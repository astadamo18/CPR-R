%CTPUTests
function result = GenVarPolyTerms(x,orders,opt)
% Construct Polynomials of Specified Orders And Part of Correction Terms
%========================================================================
% (Form)
%   result = GenVarPolyTerms(x,orders,opt)
%
% (Input)
%   x       ... Txm matrix of x_i's in columns
%   orders  ... polynomial orders (possibly cell arrays)
%   opt
%     .stoch  ... 1(default) stochastic case(calculate correction term)
%
% (Output)
%   result
%    .X     ... Tx??? matrix of x_i's polynomials
%    .Mstar ... ??? x 1 column vector of part of correction terms
%               (To make a complete Mstar, need to multiply Delta+ for each x_i
%                   where the index can be computed from result.P)
%    .P     ... (m+1)x1 vector of indicies with the first element 0, and the rest
%                elements are the ending col indicies for each x_i component in X
%=========================================================================
% <Required Files>
%   GenPowerReg.m, GenCPRCorrVec.m
%-------------------------------------------------------------------------
% (Log)
%   - Nov.15,2010: Modified for single x with nonconsecutive orders

if nargin < 3
    Stochastic = 1;
elseif isfield(opt, 'stoch')
    if opt.stoch ~= 1
        Stochastic = 0; % deterministic case
    else
        Stochastic = 1;
    end
else
    Stochastic = 1; % der default Fall, also falls kein Argument an die Funktion gegeben wird.
end
    
[T,m] = size(x); % Die Dimension der Matrix

if Stochastic == 1
    if ~iscell(orders)
        if isscalar(orders) % one max order for all x's case
            orderind = orders*ones(m,1);
        else                % diff max orders but all sequence included
            orderind = orders(:);
        end
        X = [];
        Mstar = [];
        if (m==1) && ~isscalar(orders)
            X = GenPowerReg(x,'no',orders(:));
            Mstar = GenCPRCorrVec(x,'no',orders(:));
            P = [0;length(orders)];
        else
            for i = 1:m
                X_ = GenPowerReg(x(:,i),'yes', orderind(i));
                X = [X X_];
                Mstar_ = GenCPRCorrVec(x(:,i),'yes', orderind(i));
                Mstar = [Mstar;Mstar_];
            end
            P = [0; cumsum(orderind)];   % index for the last col # for each x_i
        end
    else            % individually different orders are allowed
        X = [];
        Mstar = [];
        P = 0;
        for i = 1:m
            orderind = orders{i}(:);
            X_ = GenPowerReg(x(:,i),'no',orderind);
            X = [X X_];
            Mstar_ = GenCPRCorrVec(x(:,i),'no', orderind);
            Mstar = [Mstar;Mstar_];
            P = [P; P(end)+length(orderind)];
        end
    end
else    % Deterministic Polynomial Case
    if ~iscell(orders)
        if isscalar(orders) % one max order for all x's case
            orderind = orders*ones(m,1);
        else                % diff max orders but all sequence included
            orderind = orders(:);
        end
        X = [];
        for i = 1:m
            X_ = GenPowerReg(x(:,i),'yes', orderind(i));
            X = [X X_];
        end
        P = [0; cumsum(orderind)];   % index for the last col # for each x_i
    else            % individually different orders are allowed
        X = [];
        P = 0;
        for i = 1:m
            orderind = orders{i}(:);
            X_ = GenPowerReg(x(:,i),'no',orderind);
            X = [X X_];
            P = [P; P(end)+length(orderind)];
        end
    end
    Mstar = [];
end
result.X = X;
result.P = P;
result.Mstar = Mstar;
return

