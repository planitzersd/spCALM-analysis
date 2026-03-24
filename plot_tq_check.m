function plot_tq_check(OUT, res_tq, N, opts)
% PLOT_TQ_CHECK - Visual quality control for tq_isotonic_struct results
%
% Creates separate figures showing the raw data and estimated crossing time
% for the first N samples. Useful for quickly verifying that tq estimation
% is working as expected.
%
% INPUT
%   OUT    : structure array with .samples(i).bins (table with t, p, n)
%   res_tq : output from tq_isotonic_struct (struct array)
%   N      : number of samples to plot (default: numel(res_tq))
%   opts   : (optional) struct with fields:
%            .show_baseline  : plot horizontal lines for b and p (default true)
%            .show_threshold : plot horizontal line for theta (default true)
%            .show_windows   : shade baseline/plateau windows (default false)
%            .show_span      : overlay span_simple estimates (default false)
%            .res_span       : table from span_simple() (required if show_span=true)
%            .base_sec       : baseline window (for shading)
%            .plat_sec       : plateau window (for shading)
%            .base_frac      : baseline fraction (fallback)
%            .plat_frac      : plateau fraction (fallback)
%            .marker_size    : size of data points (default 6)
%            .save_path      : if provided, save figures to this folder
%            .close_after_save : close figures after saving (default false)
%            .save_format    : 'png'|'pdf'|'fig' (default 'png')
%
% EXAMPLE
%   % Basic usage
%   plot_tq_check(OUT, res_tq, 10);
%
%   % With window shading
%   opts = struct('show_windows', true, 'base_sec', 5, 'plat_sec', 10);
%   plot_tq_check(OUT, res_tq, 10, opts);
%
%   % Compare tq vs span estimates
%   T_span = span_simple(OUT, struct('base_sec', 5, 'plat_sec', 10));
%   opts = struct('show_span', true, 'res_span', T_span);
%   plot_tq_check(OUT, res_tq, 10, opts);
%
%   % Batch save and close (for many samples)
%   opts = struct('save_path', './qc_plots/', 'close_after_save', true);
%   plot_tq_check(OUT, res_tq, 100, opts);

%% Parse inputs
if nargin < 3 || isempty(N)
    N = min(5, numel(res_tq));
end

if nargin < 4
    opts = struct();
end

% Defaults
if ~isfield(opts, 'show_baseline'),     opts.show_baseline = true; end
if ~isfield(opts, 'show_threshold'),    opts.show_threshold = true; end
if ~isfield(opts, 'show_windows'),      opts.show_windows = false; end
if ~isfield(opts, 'show_span'),         opts.show_span = false; end
if ~isfield(opts, 'res_span'),          opts.res_span = []; end
if ~isfield(opts, 'marker_size'),       opts.marker_size = 6; end
if ~isfield(opts, 'save_path'),         opts.save_path = ''; end
if ~isfield(opts, 'close_after_save'),  opts.close_after_save = false; end
if ~isfield(opts, 'save_format'),       opts.save_format = 'png'; end
if ~isfield(opts, 'base_sec'),          opts.base_sec = []; end
if ~isfield(opts, 'plat_sec'),          opts.plat_sec = []; end
if ~isfield(opts, 'base_frac'),         opts.base_frac = 0.15; end
if ~isfield(opts, 'plat_frac'),         opts.plat_frac = 0.15; end

% Validate span option
if opts.show_span && (isempty(opts.res_span) || ~istable(opts.res_span))
    warning('show_span=true but res_span not provided or not a table. Disabling span overlay.');
    opts.show_span = false;
end

% If closing after save, must have save_path
if opts.close_after_save && isempty(opts.save_path)
    warning('close_after_save=true but no save_path provided. Figures will not be closed.');
    opts.close_after_save = false;
end

N = numel(res_tq);

%% Create save directory if needed
if ~isempty(opts.save_path) && ~exist(opts.save_path, 'dir')
    mkdir(opts.save_path);
end

%% Plot each sample
for i = 1:N
    Ss = OUT.samples(i);
    Rq = res_tq(i);
    
    % Skip if no valid bins
    if ~isfield(Ss, 'bins') || isempty(Ss.bins) || height(Ss.bins) < 1
        fprintf('Sample %d: No valid bins, skipping\n', i);
        continue;
    end
    
    Tbin = Ss.bins;
    if ~all(ismember({'t', 'p'}, Tbin.Properties.VariableNames))
        fprintf('Sample %d: Missing t or p columns, skipping\n', i);
        continue;
    end
    
    t = Tbin.t(:);
    p = Tbin.p(:);
    if ismember('n', Tbin.Properties.VariableNames)
        n = Tbin.n(:);
    else
        n = ones(size(t));
    end
    
    % Sort by time
    [t, idx] = sort(t);
    p = p(idx);
    n = n(idx);
    
    % Get span data for this sample if available
    span_row = [];
    if opts.show_span
        sample_name = get_sample_name(Ss, i);
        span_idx = find(strcmp(opts.res_span.Name, sample_name), 1);
        if ~isempty(span_idx)
            span_row = opts.res_span(span_idx, :);
        end
    end
    
    % Create new figure (invisible if closing after save)
    if opts.close_after_save
        fig = figure('Name', sprintf('Sample %d: %s', i, get_sample_name(Ss, i)), ...
                     'NumberTitle', 'off', 'Position', [100, 100, 800, 500], ...
                     'Visible', 'off');
    else
        fig = figure('Name', sprintf('Sample %d: %s', i, get_sample_name(Ss, i)), ...
                     'NumberTitle', 'off', 'Position', [100 + 50*i, 100 + 50*i, 800, 500]);
    end
    hold on; grid on;
    
    %% Shade windows if requested
    if opts.show_windows && numel(t) > 0
        % Get bin edges if available
        tlo = []; thi = [];
        if isfield(Ss, 't_lo') && isfield(Ss, 't_hi')
            if ~isempty(Ss.t_lo) && ~isempty(Ss.t_hi)
                tlo = Ss.t_lo(idx);
                thi = Ss.t_hi(idx);
            end
        end
        
        % Get window indices
        [idxB, idxP] = select_baseline_plateau_windows(t, opts, tlo, thi);
        
        % Shade baseline window
        if ~isempty(idxB)
            t_base_start = min(t(idxB));
            t_base_end = max(t(idxB));
            ylims = ylim;
            fill([t_base_start t_base_end t_base_end t_base_start], ...
                 [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                 [0.9 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        
        % Shade plateau window
        if ~isempty(idxP)
            t_plat_start = min(t(idxP));
            t_plat_end = max(t(idxP));
            ylims = ylim;
            fill([t_plat_start t_plat_end t_plat_end t_plat_start], ...
                 [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                 [1.0 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
    
    %% Plot data points
    % Size points by count if available
    if max(n) > min(n)
        scatter(t, p, opts.marker_size * (n / max(n) * 50 + 10), 'b', 'filled', ...
                'MarkerFaceAlpha', 0.6);
    else
        plot(t, p, 'b.', 'MarkerSize', opts.marker_size * 3);
    end
    
    %% Plot tq horizontal reference lines (from isotonic fit)
    if opts.show_baseline && isfinite(Rq.b)
        yline(Rq.b, '--', sprintf('tq baseline (%.3f)', Rq.b), ...
              'Color', [0.3 0.3 0.8], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    end
    
    if opts.show_baseline && isfinite(Rq.p)
        yline(Rq.p, '--', sprintf('tq plateau (%.3f)', Rq.p), ...
              'Color', [0.8 0.3 0.3], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    end
    
    if opts.show_threshold && isfinite(Rq.theta)
        yline(Rq.theta, ':', sprintf('tq threshold (%.3f)', Rq.theta), ...
              'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    end
    
    %% Plot span horizontal reference lines (from raw data)
    if opts.show_span && ~isempty(span_row)
        if isfinite(span_row.p_pre)
            yline(span_row.p_pre, '-.', sprintf('span p_{pre} (%.3f)', span_row.p_pre), ...
                  'Color', [0.5 0.5 0.9], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'right');
        end
        
        if isfinite(span_row.p_end)
            yline(span_row.p_end, '-.', sprintf('span p_{end} (%.3f)', span_row.p_end), ...
                  'Color', [0.9 0.5 0.5], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'right');
        end
    end
    
    %% Plot tq vertical line
    if isfinite(Rq.tq)
        if Rq.censored
            % Censored - dotted red line at end
            xline(Rq.tq, ':', sprintf('t_{%.0f} = %.2f (censored)', Rq.q*100, Rq.tq), ...
                  'Color', [0.8 0.2 0.2], 'LineWidth', 2, ...
                  'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        else
            % Valid crossing - solid green line
            xline(Rq.tq, '-', sprintf('t_{%.0f} = %.2f', Rq.q*100, Rq.tq), ...
                  'Color', [0.2 0.7 0.2], 'LineWidth', 2, ...
                  'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal');
        end
        
        % Add CI if available
        if isfinite(Rq.se) && Rq.se > 0
            title(sprintf('%s | t_{%.0f} = %.2f ± %.2f (SE)', ...
                  get_sample_name(Ss, i), Rq.q*100, Rq.tq, Rq.se));
        else
            title(sprintf('%s | t_{%.0f} = %.2f', ...
                  get_sample_name(Ss, i), Rq.q*100, Rq.tq));
        end
    else
        title(sprintf('%s | No valid tq estimate', get_sample_name(Ss, i)));
    end
    
    xlabel('Time');
    ylabel('Response (p)');
    
    % Add amplitude and sigma info if available
    info_str = '';
    if isfinite(Rq.A) && isfinite(Rq.sigma)
        info_str = sprintf('tq: A = %.3f, \\sigma = %.3f, A/\\sigma = %.1f', ...
                          Rq.A, Rq.sigma, Rq.A / Rq.sigma);
    end
    
    % Add span info if available
    if opts.show_span && ~isempty(span_row) && isfinite(span_row.span)
        if ~isempty(info_str)
            info_str = sprintf('%s\nspan: \\Delta = %.3f ± %.3f', ...
                              info_str, span_row.span, span_row.SE_span);
        else
            info_str = sprintf('span: \\Delta = %.3f ± %.3f', ...
                              span_row.span, span_row.SE_span);
        end
    end
    
    if ~isempty(info_str)
        text(0.02, 0.98, info_str, 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.8], ...
             'EdgeColor', 'k', 'FontSize', 9);
    end
    
    hold off;
    
    %% Save if requested
    if ~isempty(opts.save_path)
        % Determine file extension
        switch lower(opts.save_format)
            case 'png'
                ext = '.png';
            case 'pdf'
                ext = '.pdf';
            case 'fig'
                ext = '.fig';
            otherwise
                warning('Unknown save_format "%s", using png', opts.save_format);
                ext = '.png';
        end
        
        fname = fullfile(opts.save_path, sprintf('tq_check_sample_%03d%s', i, ext));
        
        try
            saveas(fig, fname);
            fprintf('Saved: %s\n', fname);
        catch ME
            warning('Failed to save figure %d: %s', i, ME.message);
        end
    end
    
    %% Close if requested
    if opts.close_after_save
        close(fig);
    end
end

fprintf('Completed plotting %d samples.\n', N);

end


%% Helper functions
function name = get_sample_name(Ss, idx)
    if isfield(Ss, 'name') && ~isempty(Ss.name)
        name = string(Ss.name);
    else
        name = sprintf('Sample %d', idx);
    end
end