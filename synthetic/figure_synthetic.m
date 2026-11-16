%% =========================================================================
% FIGURE OF THE SYNTHETIC VALIDATION (Fig. 8 of the article)
%   Top row:    metrics against the impairment profile, with the typical
%               spatial variability (path deviation of 3% of path length).
%   Bottom row: metrics against spatial variability (stress test),
%               averaged over the three profiles.
%   Sets of 5 and 8 trials; means over the replications.
%   Requires MATLAB R2020a or later (tiledlayout, exportgraphics).
% =========================================================================
close all; clear; clc;

fSynth  = fullfile('..', 'data', 'synthetic_datasets.csv');
outFile = 'figure_synthetic.pdf';
sdTyp   = 0.021;                                   % 3% of path length
sdList  = [0 0.015 0.021 0.048];                   % 0 / 2 / 3 / 6 %
sdPct   = [0 2 3 6];
profs   = ["CONTROL", "MILD", "MODERATE"];

% Methods shown (ELR with 20 landmarks is omitted: it coincides with CSR)
M = struct( ...
    'tag',   {'LIN', 'PHA', 'PVR', 'ELR3', 'CSR'}, ...
    'label', {'LIN', 'PHA', 'PVR', 'ELR ($N_R = 3$)', 'CSR'}, ...
    'color', {[0.50 0.50 0.50], [0.80 0.47 0.74], [0.90 0.60 0.00], [0.00 0.62 0.45], [0.00 0.45 0.70]}, ...
    'style', {'--', ':', '-.', '-', '-'}, ...
    'mark',  {'s', 'd', '^', 'v', 'o'});

met = {'MSD_ang_deg',     'MSD of joint angles ($^\circ$)'; ...
       'LockRet',         'Retention'; ...
       'PathErr_mean_mm', 'Path error (mm)'};

% ------------------------------- Data ------------------------------------
D = readtable(fSynth, 'TextType', 'string');
D = D(D.N_Trials >= 5, :);

% ------------------------------- Figure ----------------------------------
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 11]);
if exist('theme', 'file'), theme(fig, 'light'); end
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
letters = 'abcdef';  h = gobjects(1, numel(M));

for row = 1:2
    for k = 1:size(met, 1)
        ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
        for m = 1:numel(M)
            if row == 1                                  % against the impairment profile
                x = 1:numel(profs);  y = nan(size(x));
                for p = 1:numel(profs)
                    ix = D.Method == M(m).tag & D.Profile == profs(p) & abs(D.PathSD_m - sdTyp) < 1e-9;
                    y(p) = mean(D.(met{k, 1})(ix), 'omitnan');
                end
            else                                          % against spatial variability
                x = sdPct;  y = nan(size(x));
                for s = 1:numel(sdList)
                    ix = D.Method == M(m).tag & abs(D.PathSD_m - sdList(s)) < 1e-9;
                    y(s) = mean(D.(met{k, 1})(ix), 'omitnan');
                end
            end
            h(m) = plot(ax, x, y, M(m).style, 'Color', M(m).color, 'LineWidth', 1.3, ...
                        'Marker', M(m).mark, 'MarkerSize', 4, 'MarkerFaceColor', M(m).color);
        end
        if row == 1
            set(ax, 'XTick', 1:3, 'XTickLabel', {'Control', 'Mild', 'Moderate'}, 'XLim', [0.7 3.3]);
            xlabel(ax, 'Impairment profile', 'Interpreter', 'latex');
        else
            set(ax, 'XTick', sdPct, 'XLim', [-0.4 6.4]);
            xlabel(ax, 'Path deviation (\%)', 'Interpreter', 'latex');
        end
        ylabel(ax, met{k, 2}, 'Interpreter', 'latex');
        title(ax, sprintf('(%s)', letters(3 * (row - 1) + k)), 'Interpreter', 'latex', 'FontWeight', 'normal');
        set(ax, 'FontSize', 9, 'TickLabelInterpreter', 'latex', 'XTickLabelRotation', 0);
        yl = ylim(ax);  ylim(ax, [0 yl(2)]);
    end
end
lg = legend(h, {M.label}, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 9);
lg.Layout.Tile = 'south';

exportgraphics(fig, outFile, 'ContentType', 'vector');
fprintf('Written %s\n', outFile);
