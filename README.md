# Hand-path arc-length registration

[![DOI](https://zenodo.org/badge/1460161509.svg)](https://doi.org/10.5281/zenodo.24285469)

MATLAB code accompanying the article:

> Mayo J, García-Vallejo D, Ferrand P, Zarco MJ, Ojeda J. *Hand-path arc-length registration of post-stroke upper-limb kinematics and kinetics.* Journal of Biomechanics (submitted).

Mean curves of joint angles and moments over repeated trials require the trials to be aligned before averaging. This repository implements a registration framework that aligns trials according to the **spatial progress of the hand along its path** (the arc length of the fingertip trajectory) rather than to temporal events such as velocity peaks, together with the temporal methods used for comparison and the synthetic validation reported in the article.

## Contents

```
registration/    Functions to segment trials and normalize them with each method
  prepare_trial.m        segmentation (5% velocity threshold) and arc length of one trial
  register_trials.m      LIN, PHA, PVR, ELR and CSR normalization; mean and SD curves
  mean_sd.m              mean standard deviation (MSD) across the cycle
  peak_attenuation.m     peak attenuation of the mean curve (%)
examples/
  demo_registration.m    minimal example with simulated trials
  example_experimental.m example with the experimental marker data of one participant
synthetic/
  run_synthetic_validation.m   synthetic validation of the article (planar arm model)
  synthetic_statistics.m       statistics of the synthetic validation
  figure_synthetic.m           figure of the synthetic validation (Fig. 8)
data/
  synthetic_datasets.csv       results of the synthetic validation (one row per data set and method)
  synthetic_summary.csv        means per condition and method
  experimental/                anonymized marker trajectories of the participants (see below)
```

## Methods

| Method | Description |
|---|---|
| `LIN` | Linear time normalization from the start of the reach (IGP) to the end of the return (FRP). |
| `PHA` | Landmark registration with the four phase events (IGP, FGP, IRP, FRP). |
| `PVR` | Landmark registration with the phase events and the instants of peak velocity of the reach and the return. |
| `ELR` | Equidistant landmark registration: `NR` landmarks equally spaced in arc length in the reach and in the return. |
| `CSR` | Continuous spatial reparameterization: reach and return parameterized by normalized arc length, approach phase by normalized time, and mean curve represented against the mean time. |

In landmark registration, the time axis of each trial is warped so that its landmarks coincide with their mean across trials, using monotonic piecewise cubic Hermite interpolation. ELR with `NR = 2` is identical to PHA, and ELR converges to CSR as `NR` increases (about 20 landmarks per phase suffice in practice).

## Quick start

```matlab
addpath('registration');
cd examples
demo_registration
```

### Using your own data

For each trial you need the positions of the hand marker (e.g. the index fingertip) and the signals to be averaged (joint angles, moments, ...), sampled at the same instants:

```matlab
addpath('registration');
fs = 100;                                   % sampling frequency (Hz)
trials = cell(1, nTrials);
for k = 1:nTrials
    % pos{k}: n x 3 fingertip positions;  X{k}: n x m signals
    trials{k} = prepare_trial(pos{k}, X{k}, fs);
end
[mu, sd] = register_trials(trials, 'CSR');  % 500 x m mean and SD curves
msd = mean_sd(sd);
```

`prepare_trial` detects the reach (IGP–FGP), the approach near the target (FGP–IRP) and the return (IRP–FRP) with a threshold of 5% of the peak speed of each trial. If automatic detection fails for a trial, the events can be given manually with `opts.events = [IGP FGP IRP FRP]`. The number of output points (`K`, default 500), the internal resolution of CSR (`Np`, default 1000) and the number of ELR landmarks per phase (`NR`, default 200) can be changed with the `opts` argument of `register_trials`.

The framework assumes that the hand follows a consistent path across trials, so this should be checked before applying it to other tasks or populations.

## Reproducing the synthetic validation

```matlab
cd synthetic
run_synthetic_validation(200);      % 48 conditions x 200 replications (~45 min)
synthetic_statistics;               % reads ../data/synthetic_datasets.csv
figure_synthetic;
```

`run_synthetic_validation` writes `synthetic_datasets.csv` and `synthetic_summary.csv` in the current folder. With the default seed (1), it should reproduce the files provided in `data/`, which were generated with MATLAB. To analyse a new run, copy its output to `data/` or change the file name at the top of `synthetic_statistics.m` and `figure_synthetic.m`.

## Requirements

- **Registration functions and demo**: MATLAB R2016b or later, or GNU Octave 6 or later; no toolboxes required.
- **Synthetic validation**: Signal Processing Toolbox (`butter`, `filtfilt`), or the `signal` package in Octave.
- **Statistics and figure**: MATLAB R2020a or later; Statistics and Machine Learning Toolbox for `synthetic_statistics.m`.

## Experimental data

`data/experimental/` contains the marker trajectories of the 25 participants analysed in the article (12 controls and 13 stroke survivors: 10 mild and 3 moderate), for the trials included in the analysis. Participants are identified only by anonymous codes (C01-C12, P01-P13), the same as in Supplementary Data S1 of the article; no names, dates, clinical, demographic or anthropometric data are included.

| File | Content |
|---|---|
| `participants.csv` | Code, group (Control, Mild, Moderate), analysed side and number of trials |
| `events.csv` | Segmentation events used in the article for the analysed side (frame numbers): start (IGP) and end (FGP) of the reach, start (IRP) and end (FRP) of the return |
| `<Code>_trial<NN>.csv` | One file per trial: frame, time (s) and X, Y, Z coordinates (m) of the 16 markers |

Data were recorded at 100 Hz with an optoelectronic system (SMART-D, BTS Bioengineering) and are expressed in the laboratory reference frame. Gaps in the marker trajectories were filled and the trajectories were low-pass filtered with a fourth-order zero-lag Butterworth filter (10 Hz cut-off). Events include the manual corrections made after visual inspection. Markers: `nasion`, `rzygpr`/`lzygpr` (zygomatic processes), `rshould`/`lshould` (acromioclavicular joints), `sternum` (jugular notch), `relbow`/`lelbow` (olecranons), `rradius`/`lradius` and `rulna`/`lulna` (radial and ulnar styloids), `rhand`/`lhand` (third metacarpophalangeal joints) and `rfinger`/`lfinger` (index fingertips).

Files use commas as separators and points as decimal marks. In spreadsheet software configured with a comma decimal separator (e.g. Spanish Excel), open them with *Data > From Text/CSV* and choose the comma delimiter.

```matlab
cd examples
example_experimental('P04')        % LIN, PVR and CSR for one participant
```

The joint angles and moments of the article require the biomechanical model described in the article and are not included.

## Data availability

The marker trajectories are shared with the authorization of the Hospital Universitario Virgen del Rocío (Seville, Spain), without any information that could identify the participants. The synthetic data and all the code needed to reproduce the synthetic validation are included.

## Citation

If you use this code, please cite the article above (see `CITATION.cff`).

## License

Released under the MIT License (see `LICENSE`).

## Contact

Juana Mayo, Departamento de Ingeniería Mecánica y Fabricación, Universidad de Sevilla (juana@us.es).
