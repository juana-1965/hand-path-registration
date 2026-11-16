function demo_registration()
% DEMO_REGISTRATION  Minimal example of hand-path arc-length registration.
%
% Generates five hand-to-mouth trials of a planar two-segment arm with
% different speeds, submovements and pauses, normalizes them with linear
% time normalization (LIN), phase and peak-velocity registration (PVR),
% equidistant landmark registration (ELR) and continuous spatial
% reparameterization (CSR), and compares the inter-trial variability of the
% joint angles. Replace the generated trials with your own data: for each
% trial you need the hand marker positions and the signals to be averaged,
% sampled at the same instants.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'registration'));
rand('seed', 3);  randn('seed', 3);                        %#ok<RAND>

fs = 100;  nT = 5;
trials = cell(1, nT);
for k = 1:nT
    [pos, X] = generate_trial(fs);
    trials{k} = prepare_trial(pos, X, fs);                 % segmentation + arc length
end

methods = {'LIN', 'PVR', 'ELR', 'CSR'};
fprintf('\nInter-trial variability of the joint angles (MSD, deg)\n');
res = struct();
for i = 1:numel(methods)
    opts = struct('NR', 20);
    [mu, sd] = register_trials(trials, methods{i}, opts);
    res.(methods{i}) = struct('mu', mu, 'sd', sd);
    fprintf('  %-4s  %.2f\n', methods{i}, mean_sd(sd));
end

% Elbow angle: mean +/- SD with each method
if usejava('jvm') || exist('OCTAVE_VERSION', 'builtin')
    figure('Name', 'Hand-path arc-length registration demo');
    x = linspace(0, 100, size(res.LIN.mu, 1))';
    for i = 1:numel(methods)
        subplot(2, 2, i);  hold on;  box on;  grid on;
        r = res.(methods{i});
        fill([x; flipud(x)], [r.mu(:, 2) + r.sd(:, 2); flipud(r.mu(:, 2) - r.sd(:, 2))], ...
             [0.8 0.85 0.95], 'EdgeColor', 'none');
        plot(x, r.mu(:, 2), 'LineWidth', 1.5, 'Color', [0 0.45 0.7]);
        title(sprintf('%s (MSD %.2f deg)', methods{i}, mean_sd(r.sd)));
        xlabel('Movement cycle (%)');  ylabel('Elbow angle (deg)');
    end
end
end

% -------------------------------------------------------------------------
function [pos, X] = generate_trial(fs)
% Planar arm (upper arm 0.30 m, forearm and hand 0.35 m), shoulder at the
% origin. The fingertip goes from the table to the mouth along a curved
% path, pauses and returns. Each trial has its own durations, velocity
% profile with random submovements and pause.
A = [0.30, -0.32];  B = [0.14, 0.12];  L1 = 0.30;  L2 = 0.35;
reach = travel(A, B, 0.05, 1.2 * exp(0.15 * randn), fs);
back  = travel(B, A, 0.03, 1.1 * exp(0.15 * randn), fs);
rest  = repmat(A, round(0.3 * fs), 1);
pause = repmat(B, round((0.1 + 0.3 * rand) * fs), 1);
pos = [rest; reach; pause; back(2:end, :); rest];
pos = pos + 0.0003 * randn(size(pos));                      % marker noise
% Inverse kinematics: shoulder and elbow angles (deg)
D = (sum(pos.^2, 2) - L1^2 - L2^2) / (2 * L1 * L2);
th2 = acos(min(max(D, -1), 1));
th1 = atan2(pos(:, 2), pos(:, 1)) - atan2(L2 * sin(th2), L1 + L2 * cos(th2));
X = rad2deg([th1, th2]);
end

function p = travel(A, B, bump, T, fs)
% Curved path from A to B travelled with a bell-shaped velocity profile
% fragmented by random submovements.
n = round(T * fs) + 1;  tau = linspace(0, 1, n)';
v = tau.^2 .* (1 - tau).^2;
for j = 1:randi([0 3])                                      % submovements
    v = v .* (1 - 0.7 * exp(-0.5 * ((tau - 0.15 - 0.7 * rand) / 0.05).^2));
end
sig = cumtrapz(tau, v);  sig = sig / sig(end);
u = linspace(0, 1, 2001)';
d = (B - A) / norm(B - A);  nrm = [d(2), -d(1)];  if nrm(1) < 0, nrm = -nrm; end
xy = A + u * (B - A) + bump * sin(pi * u) * nrm;
S = [0; cumsum(sqrt(sum(diff(xy).^2, 2)))];  S = S / S(end);
p = interp1(u, xy, interp1(S, u, sig));
end
