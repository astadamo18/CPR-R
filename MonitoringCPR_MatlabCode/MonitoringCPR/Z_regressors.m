%--------------------------------------------------------------------------
% function[Z_x] = Z_regressors(Sx);
%
% Function to generate auxiliary regressors for IM-OLS regression (2)
%
% Input  Sx... Summed Regressors in R^{T \times dim(x)}
%
% Output Z_x ... Modified regressors in R^{T \times dim(x)}
%
% TJV and MW, October 2009
%--------------------------------------------------------------------------
function[Z_x] = Z_regressors(Sx);

[T,dimx] = size(Sx);

% Start now with partial summed regressors:
% Partial summed regressors:
%Sx = cumsum(x,1);

% Inverse and reinverse double sum computation:
InvBase = Sx(end:-1:1,:);
Step2 = cumsum(InvBase,1);
% Reinvert again:
RevBase = Step2(end:-1:1,:);
Z_x = cumsum(RevBase,1);

% Straightforward summation computation:
%% Time trend:
% tbase = ones(T,1);
% trend = cumsum(tbase,1);
% %
% %% Twice partial summed regressors:
% SSx = cumsum(Sx,1);
% %
% %% Triple partial summed regressors:
% SSSx = cumsum(SSx,1);
% %
% %% Overall sum of partial summed regressors:
% TotalSx = ones(1,T)*Sx;
% %
% %% Term t * Sum_1^T S_x
% T1 = trend*TotalSx;
% %%
% %
% %% Second term (double sum):
% T2 = zeros(T,dimx);
% for j=2:T,
%    T2(j,:) = SSSx(j-1,:);
% end;
% Z_x2 = T1 - T2
% %


