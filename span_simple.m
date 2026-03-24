function T = span_simple(OUT, opts)
% SPAN_SIMPLE - Simplified span estimation with optional plateau detection
%
% Estimates span = mean(p_end) - mean(p_pre) using:
%   1. Automatic plateau detection (if enabled)
%   2. Time-based windows (if base_sec/plat_sec provided)
%   3. Fraction-based windows (fallback)
%
% INPUT
%   OUT  : struct output from binned kinetics with fields:
%          OUT.samples(s).bins (table with t, p, n columns)
%          OUT.samples(s).t_lo, OUT.samples(s).t_hi (optional bin edges)
%          OUT.samples(s).name (optional)
%
%   opts : optional struct with fields:
%          % TIME-BASED windows (preferred if bin edges available)
%          .base_sec      : baseline window duration in seconds
%          .plat_sec      : plateau window duration in seconds
%          % FRACTION-BASED windows (fallback)
%          .base_frac     : fraction for early window (default 0.15)
%          .plat_frac     : fraction for late window (default 0.15)
%          % PLATEAU DETECTION
%          .detect_plateaus : attempt automatic plateau detection (default false)
%          .plateau_tol     : relative slope threshold (default 0.05)
%          .min_plateau     : min plateau size fraction (default 0.20)
%          % OTHER
%          .min_events    : min events per window (default 100)
%          .alpha         : CI level (default 0.05 for 95% CI)
%          .excel_path    : optional path for Excel export
%
% OUTPUT
%   T : table with columns Name, p_pre, p_end, span, SE_span, CIlo, CIhi,
%       method, n_pre, n_end

%% Parse options
if nargin < 2, opts = struct; end

defaults = struct('base_sec', [], ...
                  'plat_sec', [], ...
                  'base_frac', 0.05, ...
                  'plat_frac', 0.1666, ...
                  'detect_plateaus', false, ...
                  'plateau_tol', 0.05, ...
                  'min_plateau', 0.1, ...
                  'min_events', 10, ...
                  'alpha', 0.05, ...
                  'excel_path', '');

fnames = fieldnames(defaults);
for i = 1:numel(fnames)
    if ~isfield(opts, fnames{i}) || isempty(opts.(fnames{i}))
        opts.(fnames{i}) = defaults.(fnames{i});
    end
end

z = norminv(1 - opts.alpha/2);
S = numel(OUT.samples);

%% Initialize output arrays
Name = strings(S, 1);
p_pre = nan(S, 1);
p_end = nan(S, 1);
span = nan(S, 1);
SE_span = nan(S, 1);
CIlo = nan(S, 1);
CIhi = nan(S, 1);
method = strings(S, 1);
n_pre = nan(S, 1);
n_end = nan(S, 1);

%% Process each sample
for s = 1:S
    Ss = OUT.samples(s);
    
    % Get name
    if isfield(Ss, 'name') && ~isempty(Ss.name)
        Name(s) = string(Ss.name);
    else
        Name(s) = "sample_" + s;
    end
    
    % Check for valid bins
    if ~isfield(Ss, 'bins') || isempty(Ss.bins) || height(Ss.bins) < 2
        continue;
    end
    
    Tbin = Ss.bins;
    if ~all(ismember({'t', 'p', 'n'}, Tbin.Properties.VariableNames))
        continue;
    end
    
    t = Tbin.t(:);
    p = Tbin.p(:);
    n = Tbin.n(:);
    
    % Sort by time
    [t, idx] = sort(t);
    p = p(idx);
    n = n(idx);
    
    % Extract bin edges if available
    tlo = []; thi = [];
    if isfield(Ss, 't_lo') && isfield(Ss, 't_hi')
        if ~isempty(Ss.t_lo) && ~isempty(Ss.t_hi)
            tlo = Ss.t_lo(idx);
            thi = Ss.t_hi(idx);
        end
    end
    
    %% Select windows
    if opts.detect_plateaus
        % Try automatic plateau detection
        [has_early, t_early_end, has_late, t_late_start] = ...
            detect_plateaus(t, p, n, opts);
        
        if has_early && has_late
            pre_mask = t <= t_early_end;
            end_mask = t >= t_late_start;
            method(s) = "detected_plateaus";
        else
            % Fall back to window selection
            [idxB, idxP] = select_baseline_plateau_windows(t, opts, tlo, thi);
            pre_mask = false(size(t)); pre_mask(idxB) = true;
            end_mask = false(size(t)); end_mask(idxP) = true;
            method(s) = "windows_fallback";
        end
    else
        % Use shared window selection (time or fraction based)
        [idxB, idxP] = select_baseline_plateau_windows(t, opts, tlo, thi);
        pre_mask = false(size(t)); pre_mask(idxB) = true;
        end_mask = false(size(t)); end_mask(idxP) = true;
        
        if ~isempty(opts.base_sec) && ~isempty(opts.plat_sec)
            method(s) = "time_windows";
        else
            method(s) = "fraction_windows";
        end
    end
    
    %% Calculate proportions in each window
    n_pre(s) = sum(n(pre_mask));
    n_end(s) = sum(n(end_mask));
    
    % Check minimum events
    if n_pre(s) < opts.min_events || n_end(s) < opts.min_events
        continue;
    end
    
    % Weighted means (handles near-zero values correctly)
    y_pre = p(pre_mask) .* n(pre_mask);
    y_end = p(end_mask) .* n(end_mask);
    
    p_pre(s) = sum(y_pre) / n_pre(s);
    p_end(s) = sum(y_end) / n_end(s);
    
    % Standard errors (binomial)
    SE_pre = sqrt(max(p_pre(s) * (1 - p_pre(s)) / n_pre(s), 0));
    SE_end = sqrt(max(p_end(s) * (1 - p_end(s)) / n_end(s), 0));
    
    % Span and propagated error
    span(s) = p_end(s) - p_pre(s);
    SE_span(s) = sqrt(SE_pre^2 + SE_end^2);
    
    % Confidence interval
    CIlo(s) = span(s) - z * SE_span(s);
    CIhi(s) = span(s) + z * SE_span(s);
end

%% Create output table
T = table(Name, p_pre, p_end, span, SE_span, CIlo, CIhi, ...
          method, n_pre, n_end);

%% Optional Excel export
if ~isempty(opts.excel_path)
    export_to_excel(T, opts.excel_path);
end

end

%% Helper: Detect BOTH early and late plateaus
function [has_early, t_early_end, has_late, t_late_start] = detect_plateaus(t, p, n, opts)
% Detects both early (baseline) and late (response) plateaus
% Returns time boundaries for each

has_early = false;
has_late = false;
t_early_end = NaN;
t_late_start = NaN;

% Smooth the data
window = min(5, ceil(numel(p)/4));
if window < 2
    return;
end

p_smooth = movmean(p, window);

% Calculate discrete derivative
dp = diff(p_smooth);
dt = diff(t);
slope = dp ./ max(dt, eps);

% Normalize slope by max absolute value
max_slope = max(abs(slope));
if max_slope == 0
    return;
end
slope_norm = slope / max_slope;

% Find where slope drops below tolerance (plateau regions)
plateau_mask = abs(slope_norm) < opts.plateau_tol;

% Find all contiguous plateau regions
plateau_starts = find(diff([0; plateau_mask]) == 1);
plateau_ends = find(diff([plateau_mask; 0]) == -1);

if isempty(plateau_starts)
    return;
end

% Evaluate each plateau region
valid_plateaus = [];
for i = 1:numel(plateau_starts)
    pstart = plateau_starts(i);
    pend = plateau_ends(i);
    
    % Check if plateau is large enough
    plateau_frac = (pend - pstart + 1) / numel(slope);
    
    if plateau_frac >= opts.min_plateau
        % Store: [start_index, end_index, mean_time, mean_p]
        valid_plateaus = [valid_plateaus; ...
            pstart, pend, mean(t(pstart:pend)), mean(p_smooth(pstart:pend))];
    end
end

if isempty(valid_plateaus)
    return;
end

% If only one plateau found, treat it based on its position
if size(valid_plateaus, 1) == 1
    % If it's in the early part, treat as early plateau
    mid_time = (min(t) + max(t)) / 2;
    if valid_plateaus(3) < mid_time
        has_early = true;
        t_early_end = t(valid_plateaus(2));
    else
        has_late = true;
        t_late_start = t(valid_plateaus(1));
    end
    return;
end

% Multiple plateaus: pick first and last
has_early = true;
has_late = true;
t_early_end = t(valid_plateaus(1, 2));      % End of first plateau
t_late_start = t(valid_plateaus(end, 1));   % Start of last plateau

end

%% Helper: Export to Excel
function export_to_excel(T, excel_path)
xfile = string(excel_path);

if isfolder(xfile)
    xfile = fullfile(xfile, 'span_summary.xlsx');
elseif ~endsWith(xfile, {'.xlsx', '.xls'}, 'IgnoreCase', true)
    xfile = xfile + ".xlsx";
end

[folder, ~, ~] = fileparts(xfile);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end

try
    writetable(T, xfile, 'FileType', 'spreadsheet');
    fprintf('[span_simple] Exported to: %s\n', xfile);
catch ME
    warning('Excel export failed: %s', ME.message);
end
end