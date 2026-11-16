function [mu, sd, Y] = register_trials(trials, method, opts)
% REGISTER_TRIALS  Normalize a set of trials and compute their mean curve.
%
%   [mu, sd, Y] = register_trials(trials, method)
%   [mu, sd, Y] = register_trials(trials, method, opts)
%
% Inputs
%   trials  cell array of trial structures returned by prepare_trial
%   method  'LIN'  linear time normalization (IGP to FRP)
%           'PHA'  landmark registration with the four phase events
%           'PVR'  landmark registration with the phase events and the
%                  instants of peak velocity of the reach and the return
%           'ELR'  equidistant landmark registration: NR landmarks equally
%                  spaced in arc length in the reach and in the return
%           'CSR'  continuous spatial reparameterization: reach and return
%                  parameterized by normalized arc length, approach by time,
%                  mean represented against the mean time
%   opts    optional structure:
%           .K   number of points of the normalized cycle (default 500)
%           .Np  points per phase used internally by CSR (default 1000)
%           .NR  landmarks per phase for ELR (default 200)
%
% Outputs
%   mu  K x m mean curve
%   sd  K x m standard deviation across trials
%   Y   K x m x N normalized trials
%
% Landmark registration warps the time axis of each trial so that its
% landmarks coincide with their mean across trials, with monotonic
% piecewise cubic Hermite interpolation between landmarks. With NR = 2,
% ELR reduces to PHA; as NR increases, ELR converges to CSR.

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'K'),  opts.K = 500;   end
if ~isfield(opts, 'Np'), opts.Np = 1000; end
if ~isfield(opts, 'NR'), opts.NR = 200;  end

nT = numel(trials);
m = size(trials{1}.X, 2);
g = linspace(0, 1, opts.K)';
Y = zeros(opts.K, m, nT);

switch upper(method)
    case 'LIN'
        for k = 1:nT
            tr = trials{k};
            Y(:, :, k) = interp1(tr.t / tr.t(end), tr.X, g, 'pchip');
        end

    case {'PHA', 'PVR', 'ELR'}
        for k = 1:nT
            Lk = landmark_times(trials{k}, upper(method), opts.NR);
            if k == 1, L = zeros(nT, numel(Lk)); end
            L(k, :) = Lk;
        end
        Lb = mean(L, 1);
        for k = 1:nT
            tr = trials{k};
            keep = [true, diff(L(k, :)) > 0 & diff(Lb) > 0];   % drops repeated landmarks
            w = pchip(L(k, keep), Lb(keep), tr.t);              % warping function
            [w, iu] = unique(w);
            q = min(max(g * Lb(end), w(1)), w(end));
            Y(:, :, k) = interp1(w, tr.X(iu, :), q, 'pchip');
        end

    case 'CSR'
        Np = opts.Np;  q = linspace(0, 1, Np)';
        Z = zeros(3 * Np, m, nT);  Tm = zeros(3 * Np, nT);
        for k = 1:nT
            tr = trials{k};  lm = tr.lm;
            seg = {lm(1):lm(3), lm(3):lm(4), lm(4):lm(6)};    % reach, approach, return
            for p = 1:3
                ii = seg{p};  rr = (p - 1) * Np + (1:Np);
                if numel(ii) < 2                                % phase of zero duration
                    Z(rr, :, k) = repmat(tr.X(ii, :), Np, 1);
                    Tm(rr, k)   = tr.t(ii);
                    continue
                end
                if p == 2      % approach: normalized time
                    x = (tr.t(ii) - tr.t(ii(1))) / (tr.t(ii(end)) - tr.t(ii(1)));
                else           % reach and return: normalized arc length
                    x = (tr.s(ii) - tr.s(ii(1))) / (tr.s(ii(end)) - tr.s(ii(1)));
                end
                Z(rr, :, k) = interp1(x, tr.X(ii, :), q, 'pchip');
                Tm(rr, k)   = interp1(x, tr.t(ii), q, 'pchip');
            end
        end
        tb = mean(Tm, 2);  tb = tb / tb(end);                  % mean time
        [tb, iu] = unique(tb);
        for k = 1:nT
            Y(:, :, k) = interp1(tb, Z(iu, :, k), g, 'pchip');
        end

    otherwise
        error('register_trials: unknown method "%s".', method);
end

mu = mean(Y, 3);
sd = std(Y, 1, 3);
end

function L = landmark_times(tr, method, NR)
t = tr.t;  lm = tr.lm;
switch method
    case 'PHA', L = t(lm([1 3 4 6]))';
    case 'PVR', L = t(lm)';
    case 'ELR'
        r = lm(1):lm(3);  b = lm(4):lm(6);
        L = [interp1(tr.s(r), t(r), linspace(tr.s(r(1)), tr.s(r(end)), NR)), ...
             interp1(tr.s(b), t(b), linspace(tr.s(b(1)), tr.s(b(end)), NR))];
end
L = L(:)';
end
