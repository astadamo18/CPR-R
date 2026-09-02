%--------------------------------------------------------------------------
% function AMstab = AndMon_Stab(A)
% 
% Function to compute eigenvalue stabilized
% LS estimates for A(1) (AR estimate)
% Note: Discussed in paper for AR(1) pre-whitening!!!
% 
% Input:       A ... estimated AR polynomial at z=1
% Output: AMstab ... Modified version with EVs smaller equal 0.97
%--------------------------------------------------------------------------
function AMstab = AndMon_Stab(A)
[B,tmp] = eig(A*A');
[C,tmp] = eig(A'*A);
dD = diag(B'*A*C);
D = diag((abs(dD)>0.97)*0.97.*sign(dD)+(abs(dD)<=0.97).*dD);
AMstab = B*D*C';

