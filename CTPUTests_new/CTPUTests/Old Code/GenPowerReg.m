%--------------------------------------------------------------------------
% function power_reg = GenPowerReg(x,all,powvec);
%
% Function to return matrix of powers of SCALAR input regressor x.
% 
% Input:            x ... (Single) Regressor in R^{T \times 1}
%                 all ... If all == yes, then all powers from 1 to powvec 
%                         - which then is scalar - are returned.
%                         If all == no, then the powers given in powvec
%                         - then a vector - are returned.
%              powvec ... Vector in R^{1 \times NumPow}
%
% Output:   power_reg ... Matrix in R^{T \times NumPow} with powers
%--------------------------------------------------------------------------
function power_reg = GenPowerReg(x,all,powvec);

[T,sx2] = size(x);
powvec = powvec(:)';    % (2010.11.8) Added by Hong
if isequal(all,'yes');
    
    power_reg = zeros(T,powvec);
    power_reg(:,1) = x;
    for j = 2:powvec;
        power_reg(:,j) = power_reg(:,j-1).*x;
    end;
    
elseif isequal(all,'no')
    power_reg = zeros(T,size(powvec,2));
    for j = 1:size(powvec,2);
        power_reg(:,j) = x.^powvec(j);
    end;
end;


