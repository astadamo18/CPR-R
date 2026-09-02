%function result = FM_CPR(y,x,orders,w,deter,kern,band,deme,opt)

clear; clc;

%% Adatok betoltese

T_tab = readtable('panel.xlsx', 'Sheet', 'CEE');

countries = unique(T_tab.COUNTRY, 'stable');

orders = 2;
kern = 'ba';
band = 'And91';
deme = 0;
results = struct();

fprintf('%-16s %10s %10s %10s %10s %10s %10s %10s %-16s %-16s %-16s\n', ...
    'COUNTRY', 'const', 'p-value', 'GNIPC', 'p-value', 'GNIPC^2', 'p-value', "CT", "10%", "5%", "1%");

for ci = 1:numel(countries)

    cname = countries{ci};
    sub = T_tab(strcmp(T_tab.COUNTRY, cname), :);
    sub = sortrows(sub, 'YEAR');

    y = sub.NOIP / 1000;      % ezres egyseg - olvashatosag
    x = sub.GNIPC / 1000;    % ezres egyseg - numerikus stabilitas
    w = [];
    deter = ones(size(y));
    %deter = [ones(size(y)), (1:length(y))'];

    %% Sajat normcdf helyettesito (nincs szukseg Statistics Toolbox-ra)
    mynormcdf = @(z) 0.5*(1 + erf(z/sqrt(2)));
    
    FM_OLS = FM_CPR(y, x, orders, w, deter, kern, band, deme);

    const_coef  = FM_OLS.delta_fm(1);
    const_p     = 2*(1 - mynormcdf(abs(FM_OLS.t_delta(1))));

    gnipc_coef  = FM_OLS.beta_fm(1);
    gnipc_p     = 2*(1 - mynormcdf(abs(FM_OLS.t_beta(1))));

    gnipc2_coef = FM_OLS.beta_fm(2);
    gnipc2_p    = 2*(1 - mynormcdf(abs(FM_OLS.t_beta(2))));

    FM_uplus = FM_OLS.u_plus;
    FM_omega = FM_OLS.Omega_udotv1;
    d = 0;
    %d = 1;
    m = 1;
    p = 2;
    alphavec = [0.1; 0.05; 0.01];

    [CTstat, CTdecis] = CT_test(FM_uplus, FM_omega, d, m, p, alphavec);
    if CTdecis(1) == 1
    CT_10 = 'rejection';
    else
    CT_10 = 'no rejection';
    end

    if CTdecis(2) == 1
    CT_5 = 'rejection';
    else
    CT_5 = 'no rejection';
    end

    if CTdecis(3) == 1
    CT_1 = 'rejection';
    else
    CT_1 = 'no rejection';
    end

    fprintf('%-16s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %-16s %-16s %-16s\n', ...
        cname, const_coef, const_p, gnipc_coef, gnipc_p, ...
        gnipc2_coef, gnipc2_p, CTstat, CT_10, CT_5, CT_1);

    %% Eredmenyek elmentese
    fn = matlab.lang.makeValidName(cname);
    results.(fn).FM_OLS    = FM_OLS;
    results.(fn).y      = y;
    results.(fn).x      = x;
end

%% Eredmenyek mentese

save('FM_CPR_country_results.mat', 'results');