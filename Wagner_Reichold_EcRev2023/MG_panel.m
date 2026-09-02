clear; clc;
%% Adatok betoltese
T_tab = readtable('panel.xlsx', 'Sheet', 'CEE');
countries = unique(T_tab.COUNTRY, 'stable');
N = numel(countries);

%% Balanced panel epitese (kozos evek)
years = unique(T_tab.YEAR);
Y = NaN(numel(years), N);
X = NaN(numel(years), N);
for ci = 1:N
    sub = T_tab(strcmp(T_tab.COUNTRY, countries{ci}), :);
    sub = sortrows(sub, 'YEAR');
    [~, idx] = ismember(sub.YEAR, years);
    Y(idx, ci) = sub.IFDIPC / 1000;
    X(idx, ci) = sub.GNIPC / 1000;
end
if any(isnan(Y(:))) || any(isnan(X(:)))
    error('A panel nem balanced - GroupMeanFMOLS balanced panelt varva el.');
end

%% GroupMeanFMOLS futtatasa
q = 2;
type = 1;
kern = 'ba';
band = 'And91';
corrrob = 1;

[betaGM, betaiGM, V, Vdirect, alphaGM, deltaGM] = ...
    GroupMeanFMOLS(Y, X, q, type, kern, band, corrrob);

%% Orszagonkenti betak kiirasa
fprintf('%-16s %12s %12s %12s\n', 'COUNTRY','const', 'GNIPC', 'GNIPC^2');
for ci = 1:N
    fprintf('%-16s %12.3f %12.3f %12.3f\n', countries{ci}, alphaGM(ci), betaiGM(1,ci), betaiGM(2,ci));
end

fprintf('\n%-16s %12.6f %12.6f\n', 'GROUP MEAN', betaGM(1), betaGM(2));