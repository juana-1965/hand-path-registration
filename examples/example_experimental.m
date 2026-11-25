function example_experimental(code)
% EXAMPLE_EXPERIMENTAL  Apply the registration methods to the experimental
% marker data of one participant.
%
%   example_experimental('P04')
%
% Reads the trials of the participant from ../data/experimental, takes the
% fingertip marker of the analysed side as the hand path, computes a simple
% marker-based elbow flexion angle as an example signal, and compares the
% mean curves obtained with LIN, PVR and CSR. The segmentation events used
% in the article (events.csv) are passed to prepare_trial, so that the
% trials are segmented exactly as in the article.

if nargin < 1, code = 'P04'; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'registration'));
dataDir = fullfile(here, '..', 'data', 'experimental');
fs = 100;

P = read_csv_table(fullfile(dataDir, 'participants.csv'));
E = read_csv_table(fullfile(dataDir, 'events.csv'));
i = find(strcmp(P.Code, code));
if isempty(i), error('Participant %s not found in participants.csv.', code); end
s = lower(P.Analysed_side{i}(1));                      % 'r' or 'l'
nT = str2double(P.N_trials{i});

trials = {};
for k = 1:nT
    f = fullfile(dataDir, sprintf('%s_trial%02d.csv', code, k));
    if ~exist(f, 'file'), continue; end
    [D, names] = read_markers(f);
    mk = @(m) D(:, find(strcmp(names, [m '_X'])) + (0:2));
    pos = mk([s 'finger']);
    % Example signal: elbow flexion from the shoulder, elbow and wrist markers
    sh = mk([s 'should']);  el = mk([s 'elbow']);
    wr = (mk([s 'radius']) + mk([s 'ulna'])) / 2;
    a = sh - el;  b = wr - el;
    flex = 180 - acosd(sum(a .* b, 2) ./ (sqrt(sum(a.^2, 2)) .* sqrt(sum(b.^2, 2))));
    j = find(strcmp(E.Code, code) & str2double(E.Trial) == k);
    opts.events = str2double({E.IGP{j}, E.FGP{j}, E.IRP{j}, E.FRP{j}});
    trials{end + 1} = prepare_trial(pos, flex, fs, opts);     %#ok<AGROW>
end
if numel(trials) < 2, error('At least two trials are needed.'); end

methods = {'LIN', 'PVR', 'CSR'};
fprintf('\n%s (%s, %s side, %d trials): MSD of the elbow flexion angle (deg)\n', ...
        code, P.Group{i}, P.Analysed_side{i}, numel(trials));
for m = 1:numel(methods)
    [mu, sd] = register_trials(trials, methods{m});
    R.(methods{m}) = struct('mu', mu, 'sd', sd);
    fprintf('  %-4s %.2f\n', methods{m}, mean_sd(sd));
end

figure('Name', sprintf('%s: elbow flexion', code));
x = linspace(0, 100, size(R.LIN.mu, 1))';
col = {[0.5 0.5 0.5], [0.9 0.6 0], [0 0.45 0.7]};
subplot(2, 1, 1); hold on; box on; grid on;
for m = 1:numel(methods), plot(x, R.(methods{m}).mu, 'Color', col{m}, 'LineWidth', 1.4); end
ylabel('Elbow flexion (deg)'); legend(methods, 'Location', 'best'); title(code);
subplot(2, 1, 2); hold on; box on; grid on;
for m = 1:numel(methods), plot(x, R.(methods{m}).sd, 'Color', col{m}, 'LineWidth', 1.4); end
xlabel('Movement cycle (%)'); ylabel('SD across trials (deg)');
end

% -------------------------------------------------------------------------
function T = read_csv_table(f)
% Reads a small CSV file with a header into a structure of cell columns
txt = fileread(f);
lines = regexp(strtrim(txt), '\r?\n', 'split');
hdr = strsplit(strtrim(lines{1}), ',');
for c = 1:numel(hdr), T.(hdr{c}) = cell(numel(lines) - 1, 1); end
for r = 2:numel(lines)
    v = strsplit(strtrim(lines{r}), ',');
    for c = 1:numel(hdr), T.(hdr{c}){r - 1} = v{c}; end
end
end

function [D, names] = read_markers(f)
% Reads a trial file: numeric data and column names
fid = fopen(f, 'r');  names = strsplit(strtrim(fgetl(fid)), ',');  fclose(fid);
D = dlmread(f, ',', 1, 0);
end
