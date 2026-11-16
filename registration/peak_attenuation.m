function pa = peak_attenuation(mu, trials, cols)
% PEAK_ATTENUATION  Loss of peak amplitude of the mean curve (%).
%
%   pa = peak_attenuation(mu, trials, cols)
%
% For each selected column, the peak absolute value of the mean curve mu
% is compared with the average of the peak absolute values of the
% individual trials before normalization (which does not depend on the
% normalization method). The result is averaged over the columns.
if nargin < 3, cols = 1:size(mu, 2); end
pa = 0;
for c = cols
    pk = 0;
    for k = 1:numel(trials), pk = pk + max(abs(trials{k}.X(:, c))); end
    pa = pa + (1 - max(abs(mu(:, c))) / (pk / numel(trials))) * 100;
end
pa = pa / numel(cols);
end
