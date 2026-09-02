%function result = PanelEKC_indiv_eff_only(y,x,q,kern,band)

clear; clc;

%% Adatok betoltese

T_tab = readtable('panel.xlsx', 'Sheet', 'CEE');

countries = unique(T_tab.COUNTRY, 'stable');

q = 2;
kern = 'ba';
band = 'And91';
results = struct();

fprintf('%-16s %10s %10s \n', ...
    'COUNTRY', 'GNIPC', 'GNIPC^2');

for ci = 1:numel(countries)

    cname = countries{ci};
    sub = T_tab(strcmp(T_tab.COUNTRY, cname), :);
    sub = sortrows(sub, 'YEAR');

    y = sub.NOIP / 1000;      % ezres egyseg - olvashatosag
    x = sub.GNIPC / 1000;    % ezres egyseg - numerikus stabilitas
    
    FM_OLS = PanelEKC_indiv_eff_only(y, x, q, kern, band);

    gnipc_coef  = FM_OLS.beta_FM(1);
    gnipc2_coef = FM_OLS.beta_FM(2);

    fprintf('%-16s %10.4f %10.4f\n', ...
        cname, gnipc_coef, ...
        gnipc2_coef);

    %% Eredmenyek elmentese
    fn = matlab.lang.makeValidName(cname);
    results.(fn).FM_OLS    = FM_OLS;
    results.(fn).y      = y;
    results.(fn).x      = x;
end

%% Eredmenyek mentese

save('FM_CPR2_country_results.mat', 'results');