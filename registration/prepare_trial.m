function trial = prepare_trial(pos, X, fs, opts)
% PREPARE_TRIAL  Segment one hand-to-mouth trial and compute its arc length.
%
%   trial = prepare_trial(pos, X, fs)
%   trial = prepare_trial(pos, X, fs, opts)
%
% Inputs
%   pos   n x 2 or n x 3 positions of the hand (e.g. index fingertip) marker
%   X     n x m signals to be normalized (e.g. joint angles and moments),
%         sampled at the same instants as pos
%   fs    sampling frequency (Hz)
%   opts  optional structure:
%           .threshold  velocity threshold, as a fraction of the peak speed
%                       of the trial (default 0.05)
%           .minDur     minimum duration (s) for a threshold crossing to be
%                       accepted (default 0.05)
%           .events     1 x 4 sample indices [IGP FGP IRP FRP] that override
%                       the automatic segmentation (e.g. after manual review)
%
% Output: structure with the trial cropped between IGP and FRP
%   .t    time from IGP (s)
%   .X    signals
%   .s    arc length of the hand path from IGP (strictly increasing)
%   .lm   sample indices (relative to IGP) of the events
%         [IGP, peak velocity of reach, FGP, IRP, peak velocity of return, FRP]
%
% Events: IGP/FGP = start/end of the reach; IRP/FRP = start/end of the
% return. The approach phase (FGP-IRP) is the interval in which the hand
% remains near the target. If the speed does not fall below the threshold
% between both transports, FGP and IRP are set at the speed minimum in the
% central half of the movement.

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'threshold'), opts.threshold = 0.05; end
if ~isfield(opts, 'minDur'),    opts.minDur = 0.05;    end

dt = 1 / fs;
n = size(pos, 1);
vel = zeros(size(pos));
for c = 1:size(pos, 2), vel(:, c) = gradient(pos(:, c), dt); end
v = sqrt(sum(vel.^2, 2));

if isfield(opts, 'events') && ~isempty(opts.events)
    ev = opts.events(:)';
    iGP = ev(1);  iFG = ev(2);  iIR = ev(3);  iFR = ev(4);
else
    k = max(1, round(opts.minDur * fs));
    thr = opts.threshold * max(v);
    iGP = first_run(v > thr, 1, k);
    iFR = n + 1 - first_run(flipud(v > thr), 1, k);
    iFG = first_run(v < thr, iGP + 1, k);
    iIR = [];
    if ~isempty(iFG), iIR = first_run(v > thr, iFG + 1, k); end
    if isempty(iGP) || isempty(iFR) || iFR <= iGP
        error('prepare_trial: movement not detected with the current threshold.');
    end
    if isempty(iFG) || isempty(iIR) || iIR >= iFR
        % No approach phase: both transports are not separated by a pause
        a = iGP + round(0.25 * (iFR - iGP));  b = iGP + round(0.75 * (iFR - iGP));
        [~, j] = min(v(a:b));
        iFG = a + j - 1;  iIR = iFG;
    end
end

[~, j] = max(v(iGP:iFG));  iPV1 = iGP + j - 1;
[~, j] = max(v(iIR:iFR));  iPV2 = iIR + j - 1;

idx = iGP:iFR;
p = pos(idx, :);
ds = sqrt(sum(diff(p).^2, 2));
trial.t  = (0:numel(idx) - 1)' * dt;
trial.X  = X(idx, :);
trial.s  = [0; cumsum(ds)] + (0:numel(idx) - 1)' * 1e-12;   % strictly increasing
trial.lm = [iGP, iPV1, iFG, iIR, iPV2, iFR] - iGP + 1;
end

function i = first_run(b, i0, k)
% First index >= i0 from which b is true for k consecutive samples
r = conv(double(b(:)), ones(k, 1));
r = r(k:end);
i = i0 - 1 + find(r(i0:end) == k, 1, 'first');
end
