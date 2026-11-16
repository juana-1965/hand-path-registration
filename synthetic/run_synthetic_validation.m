function run_synthetic_validation(nRep, outPrefix, seed)
% RUN_SYNTHETIC_VALIDATION  Synthetic validation with known ground truth.
%
%   run_synthetic_validation(nRep, outPrefix, seed)
%     nRep      replications per condition (default 200)
%     outPrefix prefix of the output files (default 'synthetic')
%     seed      seed of the random number generator (default 1)
%
% Model: planar two-segment arm (shoulder and elbow) in the sagittal plane.
% The fingertip moves from the table to the mouth, pauses and returns. Each
% trial is generated with:
%   - its own time law (duration, timing of the velocity peak, submovements
%     that fragment the velocity profile, duration of the pause);
%   - a random spatial deformation of the hand path (stress test);
%   - a fine correction at a FIXED point of the path (a spatially
%     consistent feature present in all trials) and one at a RANDOM point;
%   - marker noise (0.3 mm), 10 Hz Butterworth filtering in both
%     directions, inverse kinematics and inverse dynamics.
%
% Ground truth: the nominal hand path, which does not depend on any time
% law. Metrics are defined in space or relative to the individual trials,
% so that no method is favoured:
%   MSD_ang_deg, MSD_mom_Nm  inter-trial variability (as in the article)
%   PathErr_mean/max_mm      distance between the hand position
%                            reconstructed from the mean joint angles
%                            (forward kinematics) and the nominal path
%   PA_ang/mom_pct           peak attenuation
%   LockRet                  retention of the fixed correction (1 = kept as
%                            in a single trial), measured by difference with
%                            paired trials without it
%
% Methods: LIN, PHA, PVR, ELR (3 and 20 landmarks per phase) and CSR.
% Phases are detected with a 5% velocity threshold that must persist for
% 50 ms.
%
% Outputs: <outPrefix>_datasets.csv (one row per data set and method) and
%          <outPrefix>_summary.csv (means per condition and method).
%
% Requires the Signal Processing Toolbox (butter, filtfilt) in MATLAB, or
% the signal package in GNU Octave.

if nargin < 1 || isempty(nRep),      nRep = 200;            end
if nargin < 2 || isempty(outPrefix), outPrefix = 'synthetic'; end
if nargin < 3 || isempty(seed),      seed = 1;              end
if exist('OCTAVE_VERSION', 'builtin'), pkg load signal; end
rng(seed);

P = parametros();
profNames  = {'CONTROL', 'MILD', 'MODERATE'};
pathSDList = [0 0.015 0.021 0.048];  % m: path deformation. Corresponds to an inter-trial
                                     % path deviation of 0 / 2 / 3 / 6 % of the path
                                     % length: none, most consistent, typical, least consistent
nTrialList = [2 3 5 8];
methNames  = {'LIN', 'PHA', 'PVR', 'ELR3', 'ELR20', 'CSR'};
metNames   = {'MSD_ang_deg', 'MSD_mom_Nm', 'PathErr_mean_mm', 'PathErr_max_mm', ...
              'PA_ang_pct', 'PA_mom_pct', 'LockRet'};
nMet = numel(metNames);

% Spatial ground truth: nominal hand path (independent of time)
[~, Ref.reachLock] = curva(P, P.A, P.B, P.bumpReach,  [0 0], P.wigLock, NaN);
[u,  Ref.reachS]   = curva(P, P.A, P.B, P.bumpReach,  [0 0], NaN,       NaN);
[~, Ref.ret]       = curva(P, P.B, P.A, P.bumpReturn, [0 0], NaN,       NaN);
Ref.u = u;
Ref.all = [Ref.reachLock; Ref.ret];

% Retention ceiling: what the nominal trial itself retains after the 10 Hz
% filter. Retention is expressed as a fraction of this ceiling (1 = ideal).
lockCeil = zeros(1, 3);
for ip = 1:3
    [e, e0] = generaEnsayo(P, P.prof.(profNames{ip}), 0, true, true);
    mN = metricas(normLin(e), zeros(500, 4), {e}, P, Ref, normLin(e0));
    lockCeil(ip) = mN(7);
end

nRows = 3 * numel(pathSDList) * numel(nTrialList) * nRep * numel(methNames);
R = nan(nRows, 5 + nMet);
row = 0;
tStart = tic;
for ip = 1:3
    pr = P.prof.(profNames{ip});
    for is = 1:numel(pathSDList)
        for in = 1:numel(nTrialList)
            nT = nTrialList(in);
            for r = 1:nRep
                tr = cell(1, nT);  tr0 = cell(1, nT);
                for k = 1:nT
                    [tr{k}, tr0{k}] = generaEnsayo(P, pr, pathSDList(is), false, true);
                end
                for m = 1:numel(methNames)
                    [mu, sd] = alinea(tr,  methNames{m}, P);
                    mu0      = alinea(tr0, methNames{m}, P);
                    row = row + 1;
                    R(row, 1:5) = [ip, is, nT, r, m];
                    R(row, 6:end) = metricas(mu, sd, tr, P, Ref, mu0);
                    R(row, end) = R(row, end) / lockCeil(ip);
                end
            end
        end
        fprintf('  %-8s pathSD = %.3f m  (%.0f s)\n', profNames{ip}, pathSDList(is), toc(tStart));
    end
end
R = R(1:row, :);

% --- One row per data set and method ---
f = fopen([outPrefix '_datasets.csv'], 'w');
fprintf(f, 'Profile,PathSD_m,N_Trials,Rep,Method');
fprintf(f, ',%s', metNames{:});  fprintf(f, '\n');
for i = 1:row
    fprintf(f, '%s,%.3f,%d,%d,%s', profNames{R(i,1)}, pathSDList(R(i,2)), R(i,3), R(i,4), methNames{R(i,5)});
    fprintf(f, ',%.6g', R(i, 6:end));  fprintf(f, '\n');
end
fclose(f);

% --- Summary per condition and method (means) ---
f = fopen([outPrefix '_summary.csv'], 'w');
fprintf(f, 'Profile,PathSD_m,N_Trials,Method,N_Datasets');
fprintf(f, ',%s', metNames{:});  fprintf(f, '\n');
for ip = 1:3
    for is = 1:numel(pathSDList)
        for in = 1:numel(nTrialList)
            for m = 1:numel(methNames)
                ix = R(:,1) == ip & R(:,2) == is & R(:,3) == nTrialList(in) & R(:,5) == m;
                fprintf(f, '%s,%.3f,%d,%s,%d', profNames{ip}, pathSDList(is), nTrialList(in), methNames{m}, sum(ix));
                fprintf(f, ',%.6g', mean(R(ix, 6:end), 1, 'omitnan'));  fprintf(f, '\n');
            end
        end
    end
end
fclose(f);
fprintf('Written %s_datasets.csv and %s_summary.csv (%.0f s)\n', outPrefix, outPrefix, toc(tStart));
end

% =========================================================================
% PARAMETERS
% =========================================================================
function P = parametros()
P.fs = 100;  P.dt = 1 / P.fs;
% Arm (segment 2 = forearm and hand, to the index fingertip)
P.L1 = 0.30;  P.L2 = 0.35;
P.m1 = 2.1;   P.m2 = 1.7;
P.lc1 = 0.13; P.lc2 = 0.17;
P.I1 = P.m1 * (0.30 * P.L1)^2;  P.I2 = P.m2 * (0.30 * P.L2)^2;
P.g = 9.81;
% Path (shoulder at the origin; x forwards, y upwards)
P.A = [0.30, -0.32];        % hand on the table
P.B = [0.14,  0.12];        % mouth
P.bumpReach  = 0.05;        % curvature of the reach (m)
P.bumpReturn = 0.03;        % curvature of the return (m)
% Fine corrections (one Hann-windowed cycle; ~5-8 Hz at the local speed)
P.wigAmp = 0.004;  P.wigWidth = 0.25;  P.wigLock = 0.85;  P.wigRandRange = [0.25 0.85];
% Measurement
P.noiseSD = 0.0003;  P.fc = 10;  [P.b, P.a] = butter(2, P.fc / (P.fs/2));
P.rest = 0.30;  P.thr = 0.05;
% Methods
P.nOut = 500;  P.Np = 1000;

% Impairment profiles, CALIBRATED against the experimental data (group means):
%   T1/T2   mean duration of the reach/return (s);  cv: inter-trial coefficient of variation
%   Tp      range of the pause duration (s)
%   pv1/pv2 mean relative timing of the velocity peak in the reach/return
%   asym    log-SD of the profile exponents (reproduces the SD of the peak timing)
%   lam     mean number of extra submovements in the reach (velocity peaks - 1)
%   lam2    idem for the return (optional; lam is used if absent)
%   depth   depth of the velocity drop of each submovement
%   wid     width of the drop (fraction of the phase)
P.prof.CONTROL  = struct('T1', 0.77, 'T2', 0.76, 'cv', 0.10, 'Tp', [0.00 0.06], 'pv1', 0.39, 'pv2', 0.44, ...
                         'asym', 0.24, 'lam', 0.02, 'depth', [0.5 0.85], 'wid', [0.04 0.08]);
P.prof.MILD     = struct('T1', 1.44, 'T2', 1.50, 'cv', 0.13, 'Tp', [0.00 0.12], 'pv1', 0.37, 'pv2', 0.20, ...
                         'asym', 0.10, 'lam', 2.30, 'lam2', 4.6, 'depth', [0.5 0.85], 'wid', [0.04 0.08]);
P.prof.MODERATE = struct('T1', 2.20, 'T2', 1.80, 'cv', 0.09, 'Tp', [0.00 0.02], 'pv1', 0.17, 'pv2', 0.22, ...
                         'asym', 0.12, 'lam', 4.30, 'lam2', 3.5, 'depth', [0.5 0.85], 'wid', [0.025 0.045]);
end

% =========================================================================
% GENERATION OF ONE TRIAL
% =========================================================================
function [tr, tr0] = generaEnsayo(P, pr, pathSD, nominal, withLock)
% tr0 (optional): paired trial with the same random draws and noise but
% without the fixed correction, used to measure retention by difference.
pq1b = 4 * [pr.pv1, 1 - pr.pv1];  pq2b = 4 * [pr.pv2, 1 - pr.pv2];   % mode = p/(p+q)
if nominal
    T1 = pr.T1;  T2 = pr.T2;  Tp = mean(pr.Tp);
    pq1 = pq1b;  pq2 = pq2b;  d1 = [];  d2 = [];  c = [0 0 0 0];  wigR = NaN;
else
    T1 = pr.T1 * exp(pr.cv * randn);  T2 = pr.T2 * exp(pr.cv * randn);
    Tp = pr.Tp(1) + diff(pr.Tp) * rand;
    pq1 = max(pq1b .* exp(pr.asym * randn(1, 2)), 1.2);   % >= 1.2: finite initial acceleration
    pq2 = max(pq2b .* exp(pr.asym * randn(1, 2)), 1.2);
    d1 = dips(pr, 1);  d2 = dips(pr, 2);
    c = pathSD * randn(1, 4);
    wigR = P.wigRandRange(1) + diff(P.wigRandRange) * rand;
end
lockPos = NaN;  if withLock, lockPos = P.wigLock; end

% Reach and return: normalized arc-length progress over time
sig1 = leyHoraria(T1, pq1, d1, P.dt);
sig2 = leyHoraria(T2, pq2, d2, P.dt);
pos2 = recorre(P, P.B, P.A, P.bumpReturn, c(3:4), NaN, NaN, sig2);
nR = round(P.rest / P.dt);  nP = max(round(Tp / P.dt), 1);
ensambla = @(p1) [repmat(P.A, nR, 1); p1; repmat(P.B, nP, 1); pos2(2:end, :); repmat(P.A, nR, 1)];

pos = ensambla(recorre(P, P.A, P.B, P.bumpReach, c(1:2), lockPos, wigR, sig1));
noise = zeros(size(pos));
if ~nominal, noise = P.noiseSD * randn(size(pos)); end
iArr = nR + numel(sig1);                  % arrival at the mouth (before cropping)
tr = procesa(P, filtfilt(P.b, P.a, pos + noise), iArr, nP);
if nargout > 1
    pos0 = ensambla(recorre(P, P.A, P.B, P.bumpReach, c(1:2), NaN, wigR, sig1));
    tr0 = procesa(P, filtfilt(P.b, P.a, pos0 + noise), iArr, nP);
end
end

function tr = procesa(P, pos, iArr, nP)

% Velocity, phases (5% threshold) and cropping to IGP-FRP
n = size(pos, 1);
vel = [gradient(pos(:, 1), P.dt), gradient(pos(:, 2), P.dt)];
v = sqrt(sum(vel.^2, 2));
thr = P.thr * max(v);
iGP = cruce(v > thr, 1);                 % crossings that persist for 50 ms (5 samples)
iFG = cruce(v < thr, iGP + 1);
iIR = cruce(v > thr, iFG + 1);
iFR = n + 1 - cruce(flipud(v > thr), 1);
iFR = min(iFR + 1, n);
% No approach phase (pause too short for velocity to fall below the threshold):
% as in the experimental data, FGP = IRP at the velocity minimum between both transports.
if isempty(iFG) || isempty(iIR) || iIR >= iFR || iFG > iArr + nP + 20
    w = max(iArr - 20, iGP + 1):min(iArr + nP + 20, iFR - 1);
    [~, a] = min(v(w));
    iFG = w(a);  iIR = iFG;
end
[~, a] = max(v(iGP:iFG));  iPV1 = iGP + a - 1;
[~, a] = max(v(iIR:iFR));  iPV2 = iIR + a - 1;

idx = iGP:iFR;
pos = pos(idx, :);
t = (0:numel(idx) - 1)' * P.dt;
lm = [iGP, iPV1, iFG, iIR, iPV2, iFR] - iGP + 1;

% Inverse kinematics and dynamics
[th1, th2] = cinInv(P, pos(:, 1), pos(:, 2));
d1t = gradient(th1, P.dt);  d2t = gradient(th2, P.dt);
dd1 = gradient(d1t, P.dt);  dd2 = gradient(d2t, P.dt);
[tau1, tau2] = dinInv(P, th1, th2, d1t, d2t, dd1, dd2);

ds = sqrt(sum(diff(pos).^2, 2));
s = [0; cumsum(ds)] + (0:numel(t) - 1)' * 1e-12;   % strictly increasing

tr.t = t;  tr.X = [rad2deg(th1), rad2deg(th2), tau1, tau2];  tr.s = s;  tr.lm = lm;
end

function i = cruce(b, i0)
% First index >= i0 from which b is true for 5 consecutive samples
k = 5;
r = conv(double(b(:)), ones(k, 1));
r = r(k:end);                             % r(i) = sum of b(i:i+k-1)
i = i0 - 1 + find(r(i0:end) == k, 1, 'first');
end

function d = dips(pr, fase)
% fase 1 = reach, 2 = return (optional lam2: submovement rate of the return)
lam = pr.lam;
if fase == 2 && isfield(pr, 'lam2'), lam = pr.lam2; end
k = poissrnd_simple(lam);
d = zeros(k, 3);
for j = 1:k
    d(j, :) = [0.12 + 0.76 * rand, pr.depth(1) + diff(pr.depth) * rand, pr.wid(1) + diff(pr.wid) * rand];
end
end

function k = poissrnd_simple(lam)
% Poisson random number without toolbox (Knuth's method)
L = exp(-lam);  k = 0;  p = rand;
while p > L, k = k + 1;  p = p * rand; end
end

function sig = leyHoraria(T, pq, d, dt)
% Beta-shaped velocity profile with multiplicative drops (submovements)
n = max(round(T / dt), 3) + 1;
tau = linspace(0, 1, n)';
v = tau.^pq(1) .* (1 - tau).^pq(2);
f = ones(n, 1);
for j = 1:size(d, 1)
    f = f .* (1 - d(j, 2) * exp(-0.5 * ((tau - d(j, 1)) / d(j, 3)).^2));
end
v = v .* max(f, 0.25);
sig = cumtrapz(tau, v);  sig = sig / sig(end);
end

function pos = recorre(P, A, B, bump, c, lockPos, wigPos, sig)
% Path A->B travelled according to the normalized arc-length progress sig.
[u, xy] = curva(P, A, B, bump, c, lockPos, wigPos);
S = [0; cumsum(sqrt(sum(diff(xy).^2, 2)))];  S = S / S(end);
uu = interp1(S, u, sig);
pos = interp1(u, xy, uu);
end

function [u, xy, nrm] = curva(P, A, B, bump, c, lockPos, wigPos)
% Path A->B with curvature, random deformation and fine corrections.
u = linspace(0, 1, 2001)';
dAB = (B - A) / norm(B - A);  nrm = [dAB(2), -dAB(1)];
if nrm(1) < 0, nrm = -nrm; end                       % forward normal
off = bump * sin(pi * u) + c(1) * sin(pi * u) + c(2) * sin(2 * pi * u);
off = off + wiggle(u, lockPos, P) + wiggle(u, wigPos, P);
xy = A + u * (B - A) + off * nrm;
end

function w = wiggle(u, u0, P)
w = zeros(size(u));
if isnan(u0), return; end
x = (u - u0) / P.wigWidth;
in = abs(x) <= 0.5;
w(in) = P.wigAmp * 0.5 * (1 + cos(2 * pi * x(in))) .* sin(2 * pi * (x(in) + 0.5));
end

function [th1, th2] = cinInv(P, x, y)
D = (x.^2 + y.^2 - P.L1^2 - P.L2^2) / (2 * P.L1 * P.L2);
th2 = acos(min(max(D, -1), 1));                      % elbow flexion (positive)
th1 = atan2(y, x) - atan2(P.L2 * sin(th2), P.L1 + P.L2 * cos(th2));
end

function [t1, t2] = dinInv(P, q1, q2, d1, d2, dd1, dd2)
c2 = cos(q2);  s2 = sin(q2);
M11 = P.I1 + P.I2 + P.m1 * P.lc1^2 + P.m2 * (P.L1^2 + P.lc2^2 + 2 * P.L1 * P.lc2 * c2);
M12 = P.I2 + P.m2 * (P.lc2^2 + P.L1 * P.lc2 * c2);
M22 = P.I2 + P.m2 * P.lc2^2;
h = P.m2 * P.L1 * P.lc2 * s2;
g1 = (P.m1 * P.lc1 + P.m2 * P.L1) * P.g * cos(q1) + P.m2 * P.lc2 * P.g * cos(q1 + q2);
g2 = P.m2 * P.lc2 * P.g * cos(q1 + q2);
t1 = M11 .* dd1 + M12 .* dd2 - h .* (2 * d1 .* d2 + d2.^2) + g1;
t2 = M12 .* dd1 + M22 .* dd2 + h .* d1.^2 + g2;
end

% =========================================================================
% NORMALIZATION METHODS
% =========================================================================
function Y = normLin(tr)
g = linspace(0, 1, 500)';
Y = interp1(tr.t / tr.t(end), tr.X, g, 'pchip');
end

function [mu, sd] = alinea(tr, meth, P)
nT = numel(tr);  g = linspace(0, 1, P.nOut)';
Y = zeros(P.nOut, 4, nT);
switch meth
    case 'LIN'
        for k = 1:nT, Y(:, :, k) = normLin(tr{k}); end
    case {'PHA', 'PVR', 'ELR3', 'ELR20'}
        for k = 1:nT
            Lk = hitos(tr{k}, meth);
            if k == 1, L = zeros(nT, numel(Lk)); end
            L(k, :) = Lk;
        end
        Lb = mean(L, 1);
        for k = 1:nT
            % Trials without approach phase (FGP = IRP): the repeated landmark
            % is omitted
            keep = [true, diff(L(k, :)) > 0 & diff(Lb) > 0];
            w = pchip(L(k, keep), Lb(keep), tr{k}.t);  % warping function
            [w, iu] = unique(w);
            qk = min(max(g * Lb(end), w(1)), w(end));   % avoids NaN from rounding at the ends
            Y(:, :, k) = interp1(w, tr{k}.X(iu, :), qk, 'pchip');
        end
    case 'CSR'
        q = linspace(0, 1, P.Np)';
        Z = zeros(3 * P.Np, 4, nT);  Tm = zeros(3 * P.Np, nT);
        for k = 1:nT
            e = tr{k};  lm = e.lm;
            seg = {lm(1):lm(3), lm(3):lm(4), lm(4):lm(6)};
            for p = 1:3
                ii = seg{p};
                rr = (p - 1) * P.Np + (1:P.Np);
                if numel(ii) < 2                      % phase of zero duration
                    Z(rr, :, k) = repmat(e.X(ii, :), P.Np, 1);
                    Tm(rr, k)   = e.t(ii);
                    continue
                end
                if p == 2
                    x = (e.t(ii) - e.t(ii(1))) / (e.t(ii(end)) - e.t(ii(1)));
                else
                    x = (e.s(ii) - e.s(ii(1))) / (e.s(ii(end)) - e.s(ii(1)));
                end
                Z(rr, :, k) = interp1(x, e.X(ii, :), q, 'pchip');
                Tm(rr, k)   = interp1(x, e.t(ii), q, 'pchip');
            end
        end
        tb = mean(Tm, 2);  tb = tb / tb(end);
        [tb, iu] = unique(tb);
        for k = 1:nT
            Y(:, :, k) = interp1(tb, Z(iu, :, k), g, 'pchip');
        end
end
mu = mean(Y, 3);
sd = std(Y, 1, 3);
end

function L = hitos(e, meth)
t = e.t;  lm = e.lm;
switch meth
    case 'PHA', L = t(lm([1 3 4 6]))';
    case 'PVR', L = t(lm)';
    otherwise
        N = str2double(meth(4:end));
        r = lm(1):lm(3);  b = lm(4):lm(6);
        L = [interp1(e.s(r), t(r), linspace(e.s(r(1)), e.s(r(end)), N)), ...
             interp1(e.s(b), t(b), linspace(e.s(b(1)), e.s(b(end)), N))];
end
L = L(:)';
end

% =========================================================================
% METRICS AGAINST THE GROUND TRUTH
% =========================================================================
function out = metricas(mu, sd, tr, P, Ref, mu0)
% Metrics that do not require defining a "true" time law:
%  - MSD (as in the article)
%  - path error: distance between the hand position reconstructed from the
%    mean joint angles (forward kinematics) and the nominal hand path.
%    Averaging misaligned postures yields postures the hand never occupied.
%  - PA, relative to the peaks of the individual trials
%  - retention of the fixed correction (1 = fully preserved)
out = nan(1, 7);
out(1) = mean(mean(sd(:, 1:2)));
out(2) = mean(mean(sd(:, 3:4)));

q1 = deg2rad(mu(:, 1));  q2 = deg2rad(mu(:, 2));
H = [P.L1 * cos(q1) + P.L2 * cos(q1 + q2), P.L1 * sin(q1) + P.L2 * sin(q1 + q2)];
dmin = zeros(size(H, 1), 1);
for j = 1:size(H, 1)
    dmin(j) = sqrt(min(sum((Ref.all - H(j, :)).^2, 2)));
end
out(3) = mean(dmin) * 1000;
out(4) = max(dmin) * 1000;

for gI = 1:2
    ch = 2 * gI - 1 + [0 1];  pa = 0;
    for c = ch
        pk = 0;
        for k = 1:numel(tr), pk = pk + max(abs(tr{k}.X(:, c))); end
        pa = pa + (1 - max(abs(mu(:, c))) / (pk / numel(tr))) * 100;
    end
    out(4 + gI) = pa / 2;
end

% Retention: difference between the mean hand path with and without the
% fixed correction (same trials, same noise), projected onto the shape of
% the correction. Deviations common to both sets cancel out.
[~, iMax] = max(mu0(:, 2));                 % end of reach ~ maximum elbow flexion
q1 = deg2rad(mu0(:, 1));  q2 = deg2rad(mu0(:, 2));
H0 = [P.L1 * cos(q1) + P.L2 * cos(q1 + q2), P.L1 * sin(q1) + P.L2 * sin(q1 + q2)];
dAB = (P.B - P.A) / norm(P.B - P.A);  nrm = [dAB(2), -dAB(1)];
if nrm(1) < 0, nrm = -nrm; end
num = 0;  den = 0;
for j = 1:iMax
    [~, k] = min(sum((Ref.reachS - H0(j, :)).^2, 2));
    wj = wiggle(Ref.u(k), P.wigLock, P);
    if wj ~= 0
        num = num + ((H(j, :) - H0(j, :)) * nrm') * wj;
        den = den + wj^2;
    end
end
if den > 0, out(7) = num / den; end
end
