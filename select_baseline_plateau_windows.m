function [idxB, idxP] = select_baseline_plateau_windows(t, opts, tlo, thi)
% SELECT_BASELINE_PLATEAU_WINDOWS
% Unified window selection for baseline and plateau regions
%
% Supports two modes:
%   1. TIME-BASED (preferred): Use actual time durations in seconds
%   2. FRACTION-BASED (fallback): Use fraction of bins/data
%
% INPUT
%   t    : [n x 1] bin centers (time values)
%   opts : struct with fields:
%          .base_sec  : baseline window duration in seconds (from start)
%          .plat_sec  : plateau window duration in seconds (to end)
%          .base_frac : baseline fraction of bins (fallback if base_sec empty)
%          .plat_frac : plateau fraction of bins (fallback if plat_sec empty)
%   tlo  : [n x 1] (optional) bin lower edges for time-based calculation
%   thi  : [n x 1] (optional) bin upper edges for time-based calculation
%
% OUTPUT
%   idxB : [1 x kB] indices for baseline bins
%   idxP : [1 x kP] indices for plateau bins
%
% EXAMPLES
%   % Time-based (if you have bin edges)
%   opts = struct('base_sec', 5, 'plat_sec', 10);
%   [idxB, idxP] = select_baseline_plateau_windows(t, opts, t_lo, t_hi);
%
%   % Fraction-based (simpler, no edges needed)
%   opts = struct('base_sec', [], 'plat_sec', [], ...
%                 'base_frac', 0.15, 'plat_frac', 0.15);
%   [idxB, idxP] = select_baseline_plateau_windows(t, opts);

    if nargin < 3, tlo = []; end
    if nargin < 4, thi = []; end
    
    % Ensure defaults
    if ~isfield(opts, 'base_sec'), opts.base_sec = []; end
    if ~isfield(opts, 'plat_sec'), opts.plat_sec = []; end
    if ~isfield(opts, 'base_frac'), opts.base_frac = 0.15; end
    if ~isfield(opts, 'plat_frac'), opts.plat_frac = 0.15; end
    
    t = t(:);
    n = numel(t);
    
    % Calculate bin widths (needed for time-based selection)
    widths = bin_widths(t, tlo, thi);
    
    % ---- BASELINE window ----
    if ~isempty(opts.base_sec)
        % Time-based: cumulative time from start
        cum = cumsum(widths);
        idxB = find(cum <= max(opts.base_sec, 0));
        if isempty(idxB), idxB = 1; end
    else
        % Fraction-based: fraction of earliest bins
        kB = max(1, ceil(opts.base_frac * n));
        idxB = 1:kB;
    end
    
    % ---- PLATEAU window ----
    if ~isempty(opts.plat_sec)
        % Time-based: cumulative time from end (going backward)
        cum_rev = cumsum(flipud(widths));
        idx_rev = find(cum_rev <= max(opts.plat_sec, 0));
        if isempty(idx_rev)
            idxP = n;
        else
            idxP = n - idx_rev + 1;
        end
    else
        % Fraction-based: fraction of latest bins
        kP = max(1, ceil(opts.plat_frac * n));
        idxP = (n - kP + 1):n;
    end
    
    % Ensure unique and row vectors
    idxB = unique(idxB(:))';
    idxP = unique(idxP(:))';
end


function widths = bin_widths(t, tlo, thi)
% BIN_WIDTHS
% Estimate per-bin time durations (in same units as t)
%
% Prefers explicit edges (tlo, thi) if available.
% Otherwise infers from bin centers using neighbor midpoints.
%
% INPUT
%   t   : [n x 1] bin centers
%   tlo : [n x 1] (optional) bin lower edges
%   thi : [n x 1] (optional) bin upper edges
%
% OUTPUT
%   widths : [n x 1] estimated duration of each bin

    n = numel(t);
    widths = zeros(n, 1);
    
    if ~isempty(tlo) && ~isempty(thi)
        % Use explicit edges
        widths = max(thi(:) - tlo(:), 0);
        
        % Fill any invalid widths with center-based estimates
        bad = ~isfinite(widths) | widths <= 0;
        if any(bad)
            widths(bad) = center_widths(t(bad), t);
        end
    else
        % Infer from centers
        widths = center_widths(t, t);
    end
    
    % Ensure positive fallback
    widths(~isfinite(widths) | widths <= 0) = eps;
end


function w = center_widths(tc, tall)
% CENTER_WIDTHS
% Compute approximate bin widths from centers using neighbor midpoints
%
% For each bin, the width is the distance between midpoints with neighbors.
% Edge bins are extrapolated symmetrically.
%
% INPUT
%   tc   : [m x 1] centers for which to compute widths
%   tall : [n x 1] all time centers (tc should be subset or equal to tall)
%
% OUTPUT
%   w : [m x 1] estimated widths

    if nargin < 2, tall = tc; end
    
    t = tall(:);
    n = numel(t);
    
    % Midpoints between consecutive bins
    mid = (t(1:end-1) + t(2:end)) / 2;
    
    % Left and right boundaries for each bin
    left = [t(1) - (mid(1) - t(1)); mid];
    right = [mid; t(end) + (t(end) - mid(end))];
    
    % Width for each bin in tall
    base_w = right - left;
    
    % Map requested centers 'tc' onto positions in 't'
    [~, loc] = ismember(tc(:), t);
    
    if any(loc == 0)
        % If tc not exactly in t, interpolate widths
        w = interp1(t, base_w, tc, 'nearest', 'extrap');
    else
        w = base_w(loc);
    end
end


function opts = parse_window_opts(opts)
% PARSE_WINDOW_OPTS
% Helper to set defaults for window selection options
%
% Ensures opts has all required fields with sensible defaults

    defaults = struct('base_sec', [], ...
                      'plat_sec', [], ...
                      'base_frac', 0.15, ...
                      'plat_frac', 0.15);
    
    fnames = fieldnames(defaults);
    for i = 1:numel(fnames)
        if ~isfield(opts, fnames{i}) || isempty(opts.(fnames{i}))
            opts.(fnames{i}) = defaults.(fnames{i});
        end
    end
end