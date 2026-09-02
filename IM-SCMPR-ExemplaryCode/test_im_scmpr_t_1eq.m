% SV, April 2024
%
% Example:
%
% Consider the following translog production function:
%
% y1 = k * c11 + l * c12 + k^2 * c13 + l^2 * c14 + (k * l) * c15 + u1
%
% with c11 = 1, c15 = 0,
% 
% i.e.:
%
% c11 = 0 0 0   c12 + 1
% c12 = 1 0 0 * c13 + 0
% c13 = 0 1 0   c14 + 0
% c14 = 0 0 1       + 0
% c15 = 0 0 0       + 0


% Generate exemplary data set:
% For sake of simplicity, here Omega = identity matrix

T = 200;
m = 2;
n = 1;

% x = cumsum(normrnd(0, 1, [T, m]));
x = cumsum(randn(T, m));
k = x(:, 1);
l = x(:, 2);

%u = normrnd(0, 1, [T, n]);
u = randn(T, n);

%cfree = normrnd(0, 1, [3, 1]);
cfree = randn(3, 1);
cfree(2) = 2;

c11 = 1;
c12 = cfree(1);
c13 = cfree(2);
c14 = cfree(3);
c15 = 0;

call = [c11 c12 c13 c14 c15];

y = [k, l, k.^2, l.^2, (k .* l)] * call' + u;


% Perform unrestricted/restricted IM estimation:

kern = 'ba';
band = 'And';
weights = eye(n);

powers = [0 1 0; 0 0 1; 0 2 0; 0 0 2; 0 1 1];

I = size(powers, 1);

D1 = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0];

D2 = eye(n * m);

d1 = [1, zeros(4, 1)']';
d2 = zeros(n * m, 1);

D = blkdiag(D1, D2);
d = [d1', d2']';

g1 = size(D1, 2); % == 3

[coeff, fitted, resid, coeff_var, resid_mod, coeff_var_mod, b, coeff_r, coeff_f, coeff_var_f] = im_scmpr(y, x, powers, D, d, weights, kern, band);


% compare coeff_f and cfree:

[coeff_f(1:g1, :), cfree]


% compare coeff_r and call:

[coeff_r(1:I, :), reshape(call', [], 1)]


% compare coeff and call:

[reshape(coeff(:, 1:I)', [], 1), reshape(call', [], 1)]


% Conventional hypothesis testing:

% Test whether c15 = 0 based on unrestricted estimates:

R = [0 0 0 0 1 0 0];
r = [0];

cent = R * reshape(coeff', [], 1) - r;

T_t = cent/sqrt(R * coeff_var * R')

%crit = norminv(0.975) # reject, if abs(T_t) > crit

% Test whether c13 = 2 based on restricted estimates:

R = [0 1 0 0 0];
r = [2];

cent = R * coeff_f - r;

T_tR = cent/sqrt(R * coeff_var_f * R')

%crit = norminv(0.975) # reject, if abs(T_tR) > crit


% Fixed-b hypothesis testing:

% Test whether c15 = 0:

R = [0 0 0 0 1 0 0];
r = [0];

cent = R * reshape(coeff', [], 1) - r;

T_tb = cent/sqrt(R * coeff_var_mod * R')

%[~, crit] = im_scmpr_fb_quantile(n, powers, size(R, 1), b, kern, 0.975) % 1 equation, 1 hypothesis per equation % reject, if abs(T_tb) > crit
