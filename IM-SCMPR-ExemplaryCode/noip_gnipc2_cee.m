% noip_gnipc2_cee.m
%
% IM-OLS becsles: NOIP = c0 + c1*GNIPC2 + c2*GNIPC2^2 + u
% orszagonkent, kulon-kulon, a panel.xlsx CEE lapjarol.
%
% FONTOS: ezt a fajlt ugyanabba a mappaba mentsd, ahol az im_scmpr.m,
% lr_var.m, lr_weights.m, And_HAC91.m, bwNW.m, commat.m talalhato,
% es a panel.xlsx is legyen elerheto (ugyanott, vagy adj meg teljes utvonalat).

clear; clc;

%% Adatok betoltese

T_tab = readtable('panel.xlsx', 'Sheet', 'CEE');

countries = unique(T_tab.COUNTRY, 'stable');

powers = [0 0; 0 1; 0 2];   % konstans, GNIPC2, GNIPC2^2
kern = 'ba';
band = 'NW';

results = struct();

fprintf('%-16s %10s %10s %10s   %10s %10s %10s\n', ...
    'Orszag', 'konst', 'GNIPC2', 'GNIPC2^2', 't_konv', 't_konv', 't_konv');

for ci = 1:numel(countries)

    cname = countries{ci};
    sub = T_tab(strcmp(T_tab.COUNTRY, cname), :);
    sub = sortrows(sub, 'YEAR');

    y = sub.NOIP / 1000;      % ezres egyseg - olvashatosag
    x = sub.GNIPC2 / 1000;    % ezres egyseg - numerikus stabilitas
    n = 1;
    m = size(x, 2);
    I = size(powers, 1);

    D = eye((I + m) * n);
    d = zeros((I + m) * n, 1);
    weights = eye(n);

    [coeff, fitted, resid, coeff_var, resid_mod, coeff_var_mod, b, ...
        coeff_r, coeff_f, coeff_var_f] = ...
        im_scmpr(y, x, powers, D, d, weights, kern, band);

    se = sqrt(diag(coeff_var));
    se_mod = sqrt(diag(coeff_var_mod));
    tstat = coeff(1:3)' ./ se(1:3);
    tstat_fb = coeff(1:3)' ./ se_mod(1:3);

    turning_point = -coeff(2) / (2 * coeff(3)) * 1000;  % nyers GNIPC2 egysegben

    fprintf('%-16s %10.4f %10.4f %10.4f   %10.2f %10.2f %10.2f\n', ...
        cname, coeff(1), coeff(2), coeff(3), tstat(1), tstat(2), tstat(3));

    results.(matlab.lang.makeValidName(cname)) = struct( ...
        'coeff', coeff, 'se', se, 'se_mod', se_mod, ...
        'tstat', tstat, 'tstat_fb', tstat_fb, ...
        'coeff_var', coeff_var, 'coeff_var_mod', coeff_var_mod, ...
        'b', b, 'T', size(y,1), 'turning_point', turning_point);

end

%% Eredmenyek mentese

save('results_cee.mat', 'results');
disp('Kesz. Eredmenyek: results struktura + results_cee.mat');