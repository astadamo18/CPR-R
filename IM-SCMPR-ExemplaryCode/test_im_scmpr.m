% SV, April 2024
%
% Example:
%
% Consider the following system of translog production functions:
%
% y1 = k * c11 + l * c12 + k^2 * c13 + l^2 * c14 + (k * l) * c15 + u1
% y2 = k * c21 + l * c22 + k^2 * c23 + l^2 * c24 + (k * l) * c25 + u2
%
% with c11 = 1, c15 = 0, c25 = 1, 
% 
% i.e., using vectorization by equations:
%
% c11 = 0 0 0 0 0 0 0   c12 + 1
% c12 = 1 0 0 0 0 0 0   c13 + 0
% c13 = 0 1 0 0 0 0 0   c14 + 0
% c14 = 0 0 1 0 0 0 0 * c21 + 0
% c15 = 0 0 0 0 0 0 0   c22 + 0
% c21 = 0 0 0 1 0 0 0   c23 + 0
% c22 = 0 0 0 0 1 0 0   c24 + 0
% c23 = 0 0 0 0 0 1 0       + 0
% c24 = 0 0 0 0 0 0 1       + 0
% c25 = 0 0 0 0 0 0 0       + 1
%
% or, using vectorization by variables:
%
% c11 = 0 0 0 0 0 0 0   c12 + 1
% c21 = 0 0 0 1 0 0 0   c13 + 0
% c12 = 1 0 0 0 0 0 0   c14 + 0
% c22 = 0 0 0 0 1 0 0 * c21 + 0
% c13 = 0 1 0 0 0 0 0   c22 + 0
% c23 = 0 0 0 0 0 1 0   c23 + 0
% c14 = 0 0 1 0 0 0 0   c24 + 0
% c24 = 0 0 0 0 0 0 1       + 0
% c15 = 0 0 0 0 0 0 0       + 0
% c25 = 0 0 0 0 0 0 0       + 1


% Generate examplary data set:
% For sake of simplicity, here Omega = identity matrix

T = 200;
m = 2;
n = 2;

x = cumsum(normrnd(0, 1, [T, m]));
k = x(:, 1);
l = x(:, 2);

u = normrnd(0, 1, [T, n]);

cfree = normrnd(0, 1, [7, 1]);
cfree(2) = 2;
cfree(4) = 3;

c11 = 1;
c12 = cfree(1);
c13 = cfree(2);
c14 = cfree(3);
c15 = 0;
c21 = cfree(4);
c22 = cfree(5);
c23 = cfree(6);
c24 = cfree(7);
c25 = 1;

call = [c11 c12 c13 c14 c15; c21 c22 c23 c24 c25];

y = [k, l, k.^2, l.^2, (k .* l)] * call' + u;


% Perform unrestricted/restricted IM estimation:

kern = 'ba';
band = 'And';
weights = eye(n);

powers = [0 1 0; 0 0 1; 0 2 0; 0 0 2; 0 1 1];

I = size(powers, 1);

D1 = [
 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0;
 1 0 0 0 0 0 0;
 0 0 0 0 1 0 0;
 0 1 0 0 0 0 0;
 0 0 0 0 0 1 0;
 0 0 1 0 0 0 0;
 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0;
];

D2 = eye(n * m);

K = commat(n, I + m);

d1 = [1, zeros(8, 1)', 1]';
d2 = zeros(n * m, 1);

D = K * blkdiag(D1, D2);
d = K * [d1', d2']';

g1 = size(D1, 2); % == 7

[coeff, fitted, resid, coeff_var, resid_mod, coeff_var_mod, b, coeff_r, coeff_f, coeff_var_f] = im_scmpr(y, x, powers, D, d, weights, kern, band);


% compare coeff_f and cfree:

[coeff_f(1:g1, :), cfree]


% compare coeff_r and call:

retain = reshape(1:((I + m) * n), [I + m, n]); 
retain = retain(1:I, :);
retain = reshape(retain, [], 1);

[coeff_r(retain, :), reshape(call', [], 1)] % == [[coeff_r(1:5, :)', coeff_r(8:12, :)']', call]

% or, using vectorization by variables:

coeff_r_variables = inv(K) * coeff_r;
call_variables = inv(commat(n, I)) * reshape(call', [], 1);

[coeff_r_variables(1:(n * I), :), call_variables]


% compare coeff and call:

% vectorization by equations:

[reshape(coeff(:, 1:I)', [], 1), reshape(call', [], 1)]

% vectorization by variables:

[reshape(coeff(:, 1:I), [], 1), reshape(call, [], 1)]


% Conventional hypothesis testing:

% Test whether c15 = 0 based on unrestricted estimates:

R = [0 0 0 0 1 0 0 0 0 0 0 0 0 0];
r = [0];

cent = R * reshape(coeff', [], 1) - r;

T_W = cent' * inv(R * coeff_var * R') * cent

%crit = chi2inv(0.95, size(R, 1))

% Test whether [c13 c21] = [2 3] based on restricted estimates:

R = [0 1 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0];
r = [2 3]';

cent = R * coeff_f - r;

T_WR = cent' * inv(R * coeff_var_f * R') * cent

%crit = chi2inv(0.95, size(R, 1))


% Fixed-b hypothesis testing:

% Test whether [c15 c25] = [0 1]:

calR = [0 0 0 0 1 0 0];
R = kron(eye(n), calR);
r = [0 1]';

cent = R * reshape(coeff', [], 1) - r;

T_Wb = cent' * inv(R * coeff_var_mod * R') * cent

%[crit, ~] = im_scmpr_fb_quantile(n, powers, size(calR, 1), b, kern, 0.95) % 2 equations, 1 hypothesis per equation

% Test whether [c15 c25] = [0 1] and [c11 c21] = [1 3]

calR = [0 0 0 0 1 0 0; 1 0 0 0 0 0 0];
R = kron(eye(n), calR);
r = [0 1 1 3]';

cent = R * reshape(coeff', [], 1) - r;

T_Wb = cent' * inv(R * coeff_var_mod * R') * cent

%[crit, ~] = im_scmpr_fb_quantile(n, powers, size(calR, 1), b, kern, 0.95) % 2 equations, 2 hypotheses per equation
