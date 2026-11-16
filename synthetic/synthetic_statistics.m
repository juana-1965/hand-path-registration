%% =========================================================================
% STATISTICS OF THE SYNTHETIC VALIDATION
% Reads ../data/synthetic_datasets.csv and writes synthetic_statistics.xlsx:
%   1_MSD_vs_PathErr  within-set association between MSD and path error
%   2_WinRates        proportion of data sets in which CSR outperforms each
%                     method, by level of spatial variability
%   3_Trials          CSR vs PVR (MSD of joint angles) by number of trials
%   4_Means           means by profile, spatial variability and method
%                     (sets of 5 and 8 trials)
% Requires the Statistics and Machine Learning Toolbox (corr).
% =========================================================================
clear; clc;
fSynth  = fullfile('..', 'data', 'synthetic_datasets.csv');
outFile = 'synthetic_statistics.xlsx';

D = readtable(fSynth, 'TextType', 'string');
meths = ["LIN", "PHA", "PVR", "ELR3", "ELR20", "CSR"];
[G, ~, Sg, Ng] = findgroups(D.Profile, D.PathSD_m, D.N_Trials, D.Rep);
[~, mIdx] = ismember(D.Method, meths);
nG = max(G);
piv = @(col) accumarray([G, mIdx], D.(col), [nG, numel(meths)], [], NaN);
iCSR = find(meths == "CSR");

% 1. Association between MSD and path error within each data set
A = piv('MSD_ang_deg');  P = piv('PathErr_mean_mm');
rho = nan(nG, 1);
for g = 1:nG
    if numel(unique(A(g, :))) > 1
        rho(g) = corr(A(g, :)', P(g, :)', 'Type', 'Spearman', 'Rows', 'complete');
    end
end
T1 = table(median(rho, 'omitnan'), mean(rho > 0, 'omitnan') * 100, ...
    'VariableNames', {'Median_rho_MSD_vs_PathErr', 'Pct_positive'});

% 2. Proportion of data sets in which CSR outperforms each method
metr = {'MSD_ang_deg', true; 'PathErr_mean_mm', true; 'LockRet', false; 'PA_mom_pct', true};
sdVals = unique(Sg);
T2 = table();
for k = 1:size(metr, 1)
    X = piv(metr{k, 1});
    for m = setdiff(1:numel(meths), iCSR)
        if metr{k, 2}, win = X(:, iCSR) < X(:, m); else, win = X(:, iCSR) > X(:, m); end
        row = table(string(metr{k, 1}), meths(m), 'VariableNames', {'Metric', 'Versus'});
        for s = 1:numel(sdVals)
            row.(sprintf('PathSD_%03dmm', round(sdVals(s) * 1000))) = mean(win(Sg == sdVals(s)));
        end
        T2 = [T2; row]; %#ok<AGROW>
    end
end

% 3. CSR vs PVR (MSD of joint angles) by number of trials
X = piv('MSD_ang_deg');  win = X(:, iCSR) < X(:, meths == "PVR");
T3 = table();
for s = 1:numel(sdVals)
    row = table(sdVals(s), 'VariableNames', {'PathSD_m'});
    for n = unique(Ng)'
        row.(sprintf('Trials_%d', n)) = mean(win(Sg == sdVals(s) & Ng == n));
    end
    T3 = [T3; row]; %#ok<AGROW>
end

% 4. Means by profile, spatial variability and method (5 and 8 trials)
Dk = D(D.N_Trials >= 5, :);
T4 = groupsummary(Dk, {'Profile', 'PathSD_m', 'Method'}, 'mean', ...
    {'MSD_ang_deg', 'MSD_mom_Nm', 'PathErr_mean_mm', 'PA_ang_pct', 'PA_mom_pct', 'LockRet'});

if isfile(outFile), delete(outFile); end
writetable(T1, outFile, 'Sheet', '1_MSD_vs_PathErr');
writetable(T2, outFile, 'Sheet', '2_WinRates');
writetable(T3, outFile, 'Sheet', '3_Trials');
writetable(T4, outFile, 'Sheet', '4_Means');
disp(T1);  fprintf('Written %s\n', outFile);
