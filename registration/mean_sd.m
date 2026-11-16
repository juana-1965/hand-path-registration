function v = mean_sd(sd, cols)
% MEAN_SD  Mean standard deviation (MSD) across the normalized cycle.
%
%   v = mean_sd(sd)        averages all columns of sd
%   v = mean_sd(sd, cols)  averages the selected columns
%
% sd is the K x m standard deviation returned by register_trials.
if nargin < 2, cols = 1:size(sd, 2); end
v = mean(mean(sd(:, cols), 1));
end
