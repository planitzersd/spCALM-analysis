function OUT = fit_binned_kinetics12(samples, opts)
%Changelog:
%   260115: added optional sample selection dialog for emission initialization
%   260112: changed mix_weight_EM calling to use p_hat from the previous
%   bin as the initialization weight for the current bin. Increased #
%   iterations. 

% ---- defaults / merge ----
if nargin<2, opts = struct; end
d = defaults();
fns = fieldnames(d);
for i=1:numel(fns)
    if ~isfield(opts, fns{i}) || isempty(opts.(fns{i})), opts.(fns{i}) = d.(fns{i}); end
end
z_alpha = norminv(1 - opts.alpha/2);

S = numel(samples);
assert(S>=1, 'No samples provided.');
for s=1:S
    samples(s).t = samples(s).t(:);
    samples(s).x = samples(s).x(:);
    assert(numel(samples(s).t)==numel(samples(s).x), 't and x mismatch for sample %d', s);
end

% ---------- 1) global emission init (with optional sample selection) ----------
if opts.emissions_use_dialog
    [params1, params2, used_idx] = auto_init_emissions_burr(samples, 'UseDialog', true, 'PlateauTimes', opts.emissions_plateau_cutoffs);
elseif ~isempty(opts.emissions_selected_samples)
    [params1, params2, used_idx] = auto_init_emissions_burr(samples, 'SelectedSamples', opts.emissions_selected_samples, 'PlateauTimes', opts.emissions_plateau_cutoffs);
else
    [params1, params2, used_idx] = auto_init_emissions_burr(samples, 'PlateauTimes', opts.emissions_plateau_cutoffs);
end

% Collect names of samples used for emission init
emissions_sample_names = cell(numel(used_idx), 1);
for i = 1:numel(used_idx)
    s = used_idx(i);
    if isfield(samples(s), 'name') && ~isempty(samples(s).name)
        emissions_sample_names{i} = samples(s).name;
    else
        emissions_sample_names{i} = sprintf('Sample %d', s);
    end
end

% ---------- 2) per-sample binning + EM ----------
OUT = struct();
OUT.params1 = params1; 
OUT.params2 = params2;
OUT.emissions_samples_used = emissions_sample_names;
OUT.emissions_indices_used = used_idx;
OUT.samples = struct([]);
OUT.opts = opts;

for s=1:S
    t = samples(s).t; x = samples(s).x;

    % per-sample target_count
    tc = choose_target_count_per_sample(numel(t), opts);

    % fill stabilizers if empty using tc
    t0_min_ct = opts.t0_edge_min_count; if isempty(t0_min_ct) || t0_min_ct==0, t0_min_ct = max(round(0.8*tc), 50); end
    t0_shr_pr = opts.t0_edge_shrink_prior; if isempty(t0_shr_pr) || t0_shr_pr<=0, t0_shr_pr = max(5, min(50, round(0.10*tc))); end

    % (optional) detect t0
    t0_est = []; se_t0 = NaN;
    if opts.t0_from_rate
        [t0_est, se_t0] = estimate_t0_from_rate( ...
            t, opts.t0_binwidth, opts.t0_smooth_w, opts.t0_direction, opts.t0_search, ...
            opts.t0_bootstrap, opts.t0_boot_reps, opts.t0_min_events_for_boot, opts.t0_alpha);
    end

    % adaptive bins
    [edges, bin] = make_bins_adaptive(t, tc, opts.max_bins, opts.min_bins, t0_est, opts.pre_t0_bins);

    % enforce minimum count near t0
    if ~isempty(t0_est) && isfinite(t0_est) && t0_min_ct > 0
        edges = enforce_min_count_near_t0(edges, t, t0_est, t0_min_ct);
        bin   = discretize(t, edges);
    end

    % per-bin EM for p
    centers = (edges(1:end-1)+edges(2:end))/2;
    nb = numel(centers);
    p_hat = nan(nb,1);
    n_i   = zeros(nb,1);
    t_lo  = edges(1:end-1).';
    t_hi  = edges(2:end).';

    for b=1:nb
        idx = (bin==b);
        xi  = x(idx);
        n_i(b) = numel(xi);
        if n_i(b)==0
            p_hat(b)=NaN;
        else
            if b == 1
                p_hat(b) = mix_weight_EM(xi, params1, params2, 0, 1000);
            else
                p_hat(b) = mix_weight_EM(xi, params1, params2, p_hat(b-1), 1000);
            end
        end
    end

    % optional shrink at t0 edge
    if ~isempty(t0_est) && isfinite(t0_est) && opts.t0_edge_shrink
        p_hat = shrink_p_near_t0(t_lo(:), t_hi(:), p_hat(:), n_i(:), t0_est, t0_shr_pr);
    end

    % ---- Robust bins table (force equal lengths & single mask)
    t_center = centers(:);
    p_col    = p_hat(:);
    n_col    = n_i(:);
    tlo_col  = t_lo(:);
    thi_col  = t_hi(:);

    keep = isfinite(t_center) & isfinite(p_col) & isfinite(n_col) & ...
           isfinite(tlo_col) & isfinite(thi_col) & (n_col > 0);

    if ~any(keep)
        bins_tbl = table( ...
            zeros(0,1), zeros(0,1), zeros(0,1,'like',n_col), zeros(0,1), zeros(0,1), ...
            'VariableNames', {'t','p','n','t_lo','t_hi'});
    else
        t_center = t_center(keep);
        p_col    = p_col(keep);
        n_col    = n_col(keep);
        tlo_col  = tlo_col(keep);
        thi_col  = thi_col(keep);

        Ls = [numel(t_center), numel(p_col), numel(n_col), numel(tlo_col), numel(thi_col)];
        assert(all(Ls == Ls(1)), 'Bins-table columns ended up with different lengths.');

        bins_tbl = table( ...
            t_center, p_col, n_col, tlo_col, thi_col, ...
            'VariableNames', {'t','p','n','t_lo','t_hi'});
    end

    OUT.samples(s).bins = bins_tbl;
    OUT.samples(s).t0   = t0_est;
    OUT.samples(s).t0_se= se_t0;

    if isfield(samples(s),'name') && ~isempty(samples(s).name)
        OUT.samples(s).name = samples(s).name;
    else
        OUT.samples(s).name = sprintf('Sample %d', s);
    end
end

% ---------- 3) fits per sample ----------
for s=1:S
    T = OUT.samples(s).bins;
    if isempty(T) || height(T)<3
        OUT.samples(s).logistic4 = [];
        OUT.samples(s).gamma     = [];
        OUT.samples(s).best      = '';
        continue
    end

    % ---- Logistic (weighted + per-sample beta selection) ----
    L = fit_bins_logistic4_optbeta( ...
            T.t, T.p, T.n, ...
            opts.display, opts.computeCI, opts.alpha, ...
            opts.logistic_weight_mode, ...
            opts.logistic_weight_beta, ...
            opts.logistic_optimize_beta, ...
            opts.logistic_beta_grid, ...
            opts.weight_floor );

    % ---- Gamma init from Logistic via quantile matching ----
    gamma_init = [];
    if ~isempty(L) && isfield(L,'params')
        B = L.params.B; A = L.params.A; t50=L.params.t50; k=L.params.k;
        t0_fixed = OUT.samples(s).t0;
        gamma_init = init_gamma_from_logistic_quintiles( ...
            T.t, B, A, t50, k, t0_fixed, opts.gamma_init_quants, opts.gamma_init_weights, opts.display);
    end

    % ---- Gamma (t0 fixed if provided; weighted + per-sample beta selection)
    G = fit_bins_gamma_onset_optbeta( ...
        T.t, T.p, T.n, ...
        opts.display, opts.computeCI, ...
        OUT.samples(s).t0, OUT.samples(s).t0_se, ...
        opts.gamma_t0_mode, ...
        opts.gamma_t0_prior_scale, ...
        opts.gamma_t0_prior_min_se, ...
        opts.gamma_t0_prior_max_se, ...
        opts.gamma_t0_prior_when_missing, ...
        opts.gamma_weight_mode, ...
        opts.gamma_weight_beta, ...
        opts.gamma_optimize_beta, ...
        opts.gamma_beta_grid, ...
        opts.weight_floor );

    OUT.samples(s).logistic4 = L;
    OUT.samples(s).gamma     = G;

    if ~isempty(L) && ~isempty(G)
        OUT.samples(s).best = ternary(L.bic <= G.bic, 'logistic4', 'gamma');
    elseif ~isempty(L)
        OUT.samples(s).best = 'logistic4';
    elseif ~isempty(G)
        OUT.samples(s).best = 'gamma';
    else
        OUT.samples(s).best = '';
    end
end

end % ===================== end main =====================


% =====================================================================
%                              DEFAULTS
% =====================================================================
function o = defaults()
o = struct();

% Binning / target count
o.target_count_mode  = 'fixed';
o.target_count       = 200;
o.target_count_scale = 0.01;
o.tc_min             = 80;

o.max_bins     = 200;
o.min_bins     = 20;
o.pre_t0_bins  = 2;

% t0 detection
o.t0_from_rate = true;
o.t0_binwidth  = 1.0;
o.t0_smooth_w  = 5;
o.t0_direction = 'decrease';
o.t0_search    = [];
o.t0_bootstrap = true;
o.t0_boot_reps = 200;
o.t0_min_events_for_boot = 200;
o.t0_alpha     = 0.05;

% Near-t0 stabilizers
o.t0_edge_min_count    = [];
o.t0_edge_shrink       = false;
o.t0_edge_shrink_prior = [];

% Weighting / CIs
o.weight_floor = 0.05;
o.alpha        = 0.05;

% Logistic weighting controls
o.logistic_weight_mode   = 'slope';
o.logistic_weight_beta   = 1;
o.logistic_optimize_beta = true;
o.logistic_beta_grid     = [0 0.5 1 1.5 2];

% Gamma weighting controls
o.gamma_weight_mode   = 'slope';
o.gamma_weight_beta   = 1;
o.gamma_optimize_beta = true;
o.gamma_beta_grid     = [0 0.5 1 1.5 2];

% Gamma t0 handling
o.gamma_t0_mode              = 'fixed';
o.gamma_t0_prior_scale       = 1.0;
o.gamma_t0_prior_min_se      = 0.05;
o.gamma_t0_prior_max_se      = Inf;
o.gamma_t0_prior_when_missing= 'free';

% Gamma init from logistic
o.gamma_init_quants  = [0.2 0.4 0.6 0.8];
o.gamma_init_weights = [];

% Emissions initialization options
o.emissions_use_dialog       = false;  % set to true to show sample selection dialog
o.emissions_selected_samples = [];     % manual indices (e.g., [1 3 5])
o.emissions_plateau_cutoffs  = [];     % user-defined time cutoffs [t_early, t_late]

% Misc
o.display      = 'off';
o.computeCI    = true;
end


% =====================================================================
%          GLOBAL EMISSIONS INIT WITH SAMPLE SELECTION
% =====================================================================
function [params1, params2, selected_idx] = auto_init_emissions_burr(samples, varargin)
% Parse optional arguments
p = inputParser;
addParameter(p, 'UseDialog', false, @islogical);
addParameter(p, 'SelectedSamples', [], @isnumeric);
addParameter(p, 'PlateauTimes', [], @isnumeric);
parse(p, varargin{:});

% Determine which samples to use
if p.Results.UseDialog
    selected_idx = show_sample_selection_dialog(samples);
    if isempty(selected_idx)
        error('No samples selected for emission initialization');
    end
else
    selected_idx = p.Results.SelectedSamples;
    if isempty(selected_idx)
        selected_idx = 1:numel(samples); % Use all samples
    end
end

% Initialize user-determined time cutoffs, if available
cutoffs = p.Results.PlateauTimes;

cap = 2e7;
early = []; late = [];
for s = selected_idx
    t = samples(s).t; x = samples(s).x;

    if isempty(cutoffs)
        t10 = quantile(t,0.10); t90 = quantile(t,0.95);
        idxE = find(t <= t10); idxL = find(t >= t90);      
    else
        if numel(cutoffs) < 2
            error('PlateauTimes must contain 2 values: [t_early, t_late]');
        end
        idxE = find(t <= cutoffs(1)); idxL = find(t >= cutoffs(2));
    end
    
    if numel(idxE) > cap, idxE = idxE(randperm(numel(idxE), cap)); end
    if numel(idxL) > cap, idxL = idxL(randperm(numel(idxL), cap)); end
    early = [early; x(idxE)]; %#ok<AGROW>
    late  = [late;  x(idxL)]; %#ok<AGROW>
    if numel(early) > 4*cap, early = early(randperm(numel(early), 2*cap)); end
    if numel(late)  > 4*cap, late  = late(randperm(numel(late),  2*cap)); end
end
assert(~isempty(early) && ~isempty(late), 'Init failed: insufficient early/late points.');

% Fit Burr distribution to early data (assumed pure State 1) - FIXED, never changes
params1 = fit_burr_distribution(early);

% Deconvolve late data if it contains mixture, then fit State 2
[late_deconv, deconv_info] = deconvolve_late_population(late, early, params1);
params2 = fit_burr_distribution(late_deconv);

% Ensure params1 corresponds to lower-valued distribution (swap if needed)
% This preserves the original params1 fit from early data
mean1 = burr_mean(params1);
mean2 = burr_mean(params2);
if mean2 < mean1
    % If State 2 mean is lower, swap the labels (shouldn't normally happen)
    warning('State 2 mean (%.2f) < State 1 mean (%.2f). Swapping labels.', mean2, mean1);
    tmp = params1; params1 = params2; params2 = tmp;
    % Also swap the data for plotting
    tmp_data = early; early = late_deconv; late_deconv = tmp_data;
end

% Plot histograms with fitted Burr distributions
plot_emission_distributions(early, late, late_deconv, params1, params2, deconv_info, selected_idx, samples);
end

function [late_deconv, info] = deconvolve_late_population(late, early, params1)
% Detect and deconvolve mixed late population using GMM and State 1 knowledge
%
% Strategy:
%   1. Check if late population appears mixed (bimodal or overlaps with early)
%   2. If mixed, fit 2-component GMM to late data
%   3. Identify which component represents converted state (State 2)
%   4. Extract high-confidence State 2 samples
%   5. If deconvolution fails or not needed, use upper percentile

late = late(:);
info = struct('method', 'none', 'fraction_kept', 1.0, 'is_mixed', false, ...
              'overlap_pct', 0, 'deconv_success', false);

% Detect mixture in late population
[is_mixed, overlap_pct] = detect_mixture_in_late(late, early);
info.is_mixed = is_mixed;
info.overlap_pct = overlap_pct;

if ~is_mixed
    % No mixture detected, use all late data
    late_deconv = late;
    info.method = 'none (pure)';
    info.fraction_kept = 1.0;
    return;
end

% Mixture detected - try GMM deconvolution
try
    % Fit 2-component GMM to late data with more iterations and better initialization
    options = statset('MaxIter', 10000, 'Display', 'off', 'TolFun', 1e-10);
    
    % Smart initialization: use k-means++ for better starting points
    gmm = fitgmdist(late, 2, 'RegularizationValue', 0.001, ...
                    'Replicates', 5, 'Options', options, 'CovarianceType', 'diagonal');
    
    % Identify which component is State 2 (converted, typically higher mean)
    % Component with higher mean is likely the converted state
    [~, state2_idx] = max(gmm.mu);
    
    % Calculate posterior probability of State 2 for each late sample
    % Correct syntax: use cluster() and convert to posterior
    cluster_idx = cluster(gmm, late);
    prob_state2 = double(cluster_idx == state2_idx);
    
    % For softer assignment, compute actual posterior probabilities
    % p(component k | x) = π_k * N(x|μ_k,Σ_k) / Σ_j π_j * N(x|μ_j,Σ_j)
    log_pdf = zeros(numel(late), 2);
    for k = 1:2
        log_pdf(:, k) = log(gmm.ComponentProportion(k)) + ...
                        log(normpdf(late, gmm.mu(k), sqrt(gmm.Sigma(k))));
    end
    % Numerically stable softmax
    log_pdf_max = max(log_pdf, [], 2);
    pdf_norm = exp(log_pdf - log_pdf_max);
    posterior_probs = pdf_norm ./ sum(pdf_norm, 2);
    prob_state2 = posterior_probs(:, state2_idx);
    
    % Extract high-confidence State 2 samples
    confidence_threshold = 0.35;  % Lowered from 0.6 for heavily mixed data
    high_conf_idx = prob_state2 > confidence_threshold;
    
    % Require minimum samples for reliable fit
    min_samples = max(100, 0.05 * numel(late));  % Lowered from 0.1
    
    if sum(high_conf_idx) >= min_samples
        late_deconv = late(high_conf_idx);
        info.method = 'GMM deconvolution';
        info.fraction_kept = sum(high_conf_idx) / numel(late);
        info.deconv_success = true;
        info.gmm_mixing = gmm.ComponentProportion(state2_idx);
        info.confidence_threshold = confidence_threshold;
    else
        % Not enough high-confidence samples, try weighted approach
        late_deconv = weighted_percentile_deconv(late, early, prob_state2);
        info.method = 'GMM-weighted percentile';
        info.fraction_kept = numel(late_deconv) / numel(late);
    end
    
catch ME
    % GMM failed, fall back to percentile method
    warning('GMM deconvolution failed: %s. Using percentile method.', ME.message);
    late_deconv = percentile_based_deconv(late, early);
    info.method = 'percentile (GMM failed)';
    info.fraction_kept = numel(late_deconv) / numel(late);
end

% Ensure we have enough data
if numel(late_deconv) < 50
    warning('Deconvolution resulted in too few samples (%d). Using upper 40%% of late data.', ...
            numel(late_deconv));
    late_deconv = late(late > quantile(late, 0.60));
    info.method = 'percentile fallback';
    info.fraction_kept = 0.4;
end
end

function [is_mixed, overlap_pct] = detect_mixture_in_late(late, early)
% Detect if late population is mixed with unconverted state
%
% Criteria:
%   1. Significant overlap with early distribution
%   2. Bimodality indicators

late = late(:);
early = early(:);

% Criterion 1: Overlap with early distribution
early_upper = quantile(early, 0.90);  % upper bound of early state
overlap_count = sum(late <= early_upper);
overlap_pct = 100 * overlap_count / numel(late);

% Significant overlap if >20% of late data falls within early range
high_overlap = overlap_pct > 20;

% Criterion 2: Bimodality check using simple second derivative of histogram
% (lighter weight than Hartigan's dip test)
edges = linspace(min(late), max(late), 50);
counts = histcounts(late, edges);
counts_smooth = movmean(counts, 3);  % smooth to reduce noise

% Find local maxima (peaks)
d2 = diff(diff([0 counts_smooth 0]));  % second derivative
peaks = find(d2(1:end-1) < 0 & d2(2:end) > 0);
has_multiple_peaks = numel(peaks) > 1;

% Consider mixed if either criterion met
is_mixed = high_overlap || has_multiple_peaks;
end

function late_deconv = percentile_based_deconv(late, early)
% Fallback: use upper percentile of late data that's beyond early range

late = late(:);

% Method 1: Use data beyond early distribution's upper tail
early_upper = quantile(early, 0.90);  % Lowered from 0.95 for better separation
beyond_early = late(late > early_upper);

if numel(beyond_early) >= 100
    late_deconv = beyond_early;
else
    % Method 2: Use upper portion of late data
    % For heavily mixed data, be more aggressive
    late_threshold = quantile(late, 0.50);  % Use upper 50%
    late_deconv = late(late > late_threshold);
end

% Ensure minimum sample size
if numel(late_deconv) < 100 && numel(late) >= 100
    late_threshold = quantile(late, max(0.4, 1 - 100/numel(late)));
    late_deconv = late(late > late_threshold);
end
end

function late_deconv = weighted_percentile_deconv(late, early, weights)
% Use GMM weights to extract most likely State 2 samples
% Take top samples by GMM probability

late = late(:);
weights = weights(:);

% Sort by probability of being State 2
[sorted_weights, sort_idx] = sort(weights, 'descend');

% Take top samples representing at least 15% of data or until weights drop below 0.3
cumsum_weights = cumsum(sorted_weights) / sum(sorted_weights);
cutoff_by_mass = find(cumsum_weights >= 0.15, 1, 'first');
cutoff_by_prob = find(sorted_weights < 0.3, 1, 'first') - 1;

if isempty(cutoff_by_prob), cutoff_by_prob = numel(late); end
if isempty(cutoff_by_mass), cutoff_by_mass = round(0.15 * numel(late)); end

n_keep = max(min(cutoff_by_mass, cutoff_by_prob), 100);  % At least 100 samples
n_keep = min(n_keep, numel(late));

late_deconv = late(sort_idx(1:n_keep));
end

function plot_emission_distributions(early, late, late_deconv, params1, params2, deconv_info, selected_idx, samples)
% Create figure with overlaid histograms and fitted Burr distributions
figure('Name', 'Emission Distribution Initialization', 'Position', [100 100 800 500]);

% Determine shared x-axis range
all_data = [early(:); late(:)];
x_min = min(all_data);
x_max = max(all_data);
x_range = x_max - x_min;
x_plot = linspace(max(x_min - 0.1*x_range, eps), x_max + 0.1*x_range, 500);

% Plot overlaid histograms
h1 = histogram(early, 50, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
h2 = histogram(late, 50, 'Normalization', 'pdf', 'FaceColor', [0.9 0.5 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% If deconvolution was performed, show the deconvolved data
if deconv_info.is_mixed && ~strcmp(deconv_info.method, 'none (pure)')
    h3 = histogram(late_deconv, 40, 'Normalization', 'pdf', 'FaceColor', [0.8 0.2 0.2], ...
                   'EdgeColor', [0.5 0 0], 'LineWidth', 1.5, 'FaceAlpha', 0.6);
end

% Plot fitted Burr distributions
pdf1 = pdf('Burr', x_plot, params1(1), params1(2), params1(3));
pdf2 = pdf('Burr', x_plot, params2(1), params2(2), params2(3));
plot(x_plot, pdf1, 'b-', 'LineWidth', 2.5);
plot(x_plot, pdf2, 'r-', 'LineWidth', 2.5);

% NEW: Plot mixture model reconstruction if deconvolution occurred
if deconv_info.is_mixed && ~strcmp(deconv_info.method, 'none (pure)')
    % Estimate mixing proportion from deconvolution
    p_state2 = deconv_info.fraction_kept;  % fraction identified as State 2
    p_state1 = 1 - p_state2;                % remaining is State 1
    
    % Mixture PDF: p_state1 * pdf1 + p_state2 * pdf2
    pdf_mixture = p_state1 * pdf1 + p_state2 * pdf2;
    plot(x_plot, pdf_mixture, 'm-', 'LineWidth', 2.5, 'LineStyle', '--');
    
    % Calculate reconstruction error (KL divergence or simple L2)
    % Compare mixture to actual late data histogram
    [late_counts, late_edges] = histcounts(late, 50, 'Normalization', 'pdf');
    late_centers = (late_edges(1:end-1) + late_edges(2:end)) / 2;
    
    % Interpolate mixture PDF at histogram centers
    pdf_mixture_interp = interp1(x_plot, pdf_mixture, late_centers, 'linear', 0);
    
    % L2 norm (RMSE)
    valid_idx = isfinite(late_counts) & isfinite(pdf_mixture_interp);
    if sum(valid_idx) > 0
        rmse = sqrt(mean((late_counts(valid_idx) - pdf_mixture_interp(valid_idx)).^2));
        mixture_quality = 1 - min(rmse / mean(late_counts(valid_idx)), 1);  % 0 to 1 scale
    else
        mixture_quality = NaN;
    end
else
    mixture_quality = NaN;
    p_state1 = NaN;
    p_state2 = NaN;
end

xlabel('Fluorescence', 'FontSize', 11);
ylabel('Probability Density', 'FontSize', 11);
title('Emission Distribution Initialization', 'FontSize', 12, 'FontWeight', 'bold');

% Create legend with parameter info and deconvolution status
if deconv_info.is_mixed && ~strcmp(deconv_info.method, 'none (pure)')
    legend({'Early Data (State 1)', 'Late Data (raw)', 'Late Data (deconvolved)', ...
            sprintf('State 1: Burr(α=%.3f, c=%.2f, k=%.3f)', params1(1), params1(2), params1(3)), ...
            sprintf('State 2: Burr(α=%.3f, c=%.2f, k=%.3f)', params2(1), params2(2), params2(3)), ...
            sprintf('Mixture: %.0f%% S1 + %.0f%% S2', 100*p_state1, 100*p_state2)}, ...
           'Location', 'best', 'FontSize', 9);
else
    legend({'Early Data (State 1)', 'Late Data (State 2)', ...
            sprintf('State 1: Burr(α=%.3f, c=%.2f, k=%.3f)', params1(1), params1(2), params1(3)), ...
            sprintf('State 2: Burr(α=%.3f, c=%.2f, k=%.3f)', params2(1), params2(2), params2(3))}, ...
           'Location', 'best', 'FontSize', 9);
end

grid on;
set(gca, 'FontSize', 10);

% Add sample information and deconvolution details
sample_names_str = get_sample_names_string(selected_idx, samples);
info_lines = {sprintf('Samples: %s', sample_names_str)};

if deconv_info.is_mixed
    info_lines{end+1} = sprintf('Mixture detected (%.0f%% overlap)', deconv_info.overlap_pct);
    info_lines{end+1} = sprintf('Method: %s', deconv_info.method);
    info_lines{end+1} = sprintf('Kept: %.0f%% of late data', 100*deconv_info.fraction_kept);
    if isfinite(mixture_quality)
        info_lines{end+1} = sprintf('Reconstruction quality: %.1f%%', 100*mixture_quality);
    end
end

info_text = strjoin(info_lines, '\n');
text(0.02, 0.98, info_text, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 9, 'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.7 0.7 0.7]);

drawnow;
end

function str = get_sample_names_string(selected_idx, samples)
% Create a compact string of sample names for the title
if numel(selected_idx) <= 5
    % Show all names if 5 or fewer
    names = cell(1, numel(selected_idx));
    for i = 1:numel(selected_idx)
        s = selected_idx(i);
        if isfield(samples(s), 'name') && ~isempty(samples(s).name)
            names{i} = samples(s).name;
        else
            names{i} = sprintf('S%d', s);
        end
    end
    str = strjoin(names, ', ');
else
    % Show first few and count if many
    str = sprintf('%d samples', numel(selected_idx));
end
end

function selected_idx = show_sample_selection_dialog(samples)
% Create sample names list
sample_names = cell(numel(samples), 1);
for i = 1:numel(samples)
    if isfield(samples(i), 'name') && ~isempty(samples(i).name)
        sample_names{i} = samples(i).name;
    else
        sample_names{i} = sprintf('Sample %d', i);
    end
end

% Create dialog figure
fig = dialog('Name', 'Select Samples for Emission Initialization', ...
             'Position', [100 100 450 550]);

% Add instruction text
uicontrol('Parent', fig, 'Style', 'text', ...
          'String', 'Select samples to use for initializing emission parameters:', ...
          'Position', [20 510 410 30], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left');

% Add help text
uicontrol('Parent', fig, 'Style', 'text', ...
          'String', 'Choose samples representing clear early/late states', ...
          'Position', [20 480 410 20], ...
          'FontSize', 9, ...
          'ForegroundColor', [0.4 0.4 0.4], ...
          'HorizontalAlignment', 'left');

% Create listbox with multiple selection
listbox = uicontrol('Parent', fig, 'Style', 'listbox', ...
                    'String', sample_names, ...
                    'Min', 0, 'Max', numel(samples), ...
                    'Position', [20 130 410 340], ...
                    'Value', 1:numel(samples)); % All selected by default

% Add select all button
uicontrol('Parent', fig, 'Style', 'pushbutton', ...
          'String', 'Select All', ...
          'Position', [20 85 100 30], ...
          'Callback', @(~,~) set(listbox, 'Value', 1:numel(samples)));

% Add clear all button
uicontrol('Parent', fig, 'Style', 'pushbutton', ...
          'String', 'Clear All', ...
          'Position', [130 85 100 30], ...
          'Callback', @(~,~) set(listbox, 'Value', []));

% Status text
status_txt = uicontrol('Parent', fig, 'Style', 'text', ...
          'String', sprintf('%d of %d samples selected', numel(samples), numel(samples)), ...
          'Position', [20 55 410 20], ...
          'FontSize', 9, ...
          'HorizontalAlignment', 'center');

% Update status on selection change
set(listbox, 'Callback', @(src,~) set(status_txt, 'String', ...
    sprintf('%d of %d samples selected', numel(get(src,'Value')), numel(samples))));

% Add OK button
selected_idx = [];
uicontrol('Parent', fig, 'Style', 'pushbutton', ...
          'String', 'OK', ...
          'Position', [250 15 90 30], ...
          'Callback', @ok_callback);

% Add Cancel button
uicontrol('Parent', fig, 'Style', 'pushbutton', ...
          'String', 'Cancel', ...
          'Position', [350 15 80 30], ...
          'Callback', @cancel_callback);

% Wait for user input
uiwait(fig);

    function ok_callback(~, ~)
        selected_idx = get(listbox, 'Value');
        if isempty(selected_idx)
            warndlg('Please select at least one sample', 'No Selection');
        else
            delete(fig);
        end
    end

    function cancel_callback(~, ~)
        selected_idx = [];
        delete(fig);
    end
end

function params = fit_burr_distribution(data)
data = data(data > 0);
assert(~isempty(data), 'No positive data points for Burr fitting');

% Try multiple approaches in order of preference
params = [];

% Approach 1: MLE with relaxed constraints
try
    % Add lower bounds to prevent degenerate solutions
    % Burr parameters: [alpha, c, k] where alpha > 0, c > 0, k > 0
    options = statset('MaxIter', 2000, 'MaxFunEvals', 5000, 'Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8);
    
    % Starting values based on data
    alpha_start = median(data);
    c_start = 3.0;
    k_start = 2.0;
    
    % Use custom MLE with bounds (relaxed upper bounds)
    params_fit = mle(data, 'pdf', @(x, a, c, k) pdf('Burr', x, a, c, k), ...
                     'start', [alpha_start, c_start, k_start], ...
                     'lower', [eps, 0.1, 0.01], ...
                     'upper', [max(data)*10, 500, 100], ...
                     'options', options);
    
    % Validate the fit - reject if still hitting bounds or extreme
    if all(isfinite(params_fit)) && all(params_fit > 0) && ...
       params_fit(2) < 450 && params_fit(3) < 90 && params_fit(3) > 0.02
        params = params_fit;
        return;
    end
catch ME
    % MLE failed, try next approach
end

% Approach 2: Standard MLE without constraints
if isempty(params)
    try
        options2 = statset('MaxIter', 2000, 'Display', 'off');
        params_fit = mle(data, 'distribution', 'Burr', 'Options', options2);
        
        % Validate: reject if parameters are extreme
        if all(isfinite(params_fit)) && all(params_fit > 0) && ...
           params_fit(2) < 1000 && params_fit(3) < 1000 && params_fit(3) > 0.001
            params = params_fit;
            return;
        end
    catch
        % Continue to fallback
    end
end

% Approach 3: Method of moments / percentile matching
if isempty(params)
    warning('Burr MLE failed, using robust moment-based estimate');
    
    % Use robust statistics
    alpha = median(data);
    
    % Estimate c and k from shape
    q10 = quantile(data, 0.10);
    q25 = quantile(data, 0.25);
    q50 = median(data);
    q75 = quantile(data, 0.75);
    q90 = quantile(data, 0.90);
    
    % IQR-based shape estimation
    iqr_val = q75 - q25;
    skew_proxy = (q75 - q50) / max(q50 - q25, eps);  % right skew indicator
    
    % Kurtosis proxy from tail behavior
    tail_ratio = (q90 - q50) / max(q50 - q10, eps);
    
    % Set c based on spread (lower c = heavier tails)
    c = max(1.0, min(20, 5 / max(iqr_val / q50, 0.1)));  % reasonable range
    
    % Set k based on skewness (lower k = more right-skewed)
    k = max(0.5, min(10, 2.5 / max(skew_proxy, 0.3)));
    
    params = [alpha, c, k];
end

% Final validation
assert(all(isfinite(params)) && all(params > 0), 'Burr parameter estimation failed');
end

function m = burr_mean(params)
alpha = params(1);
c = params(2);
k = params(3);

if k > 1/c
    m = alpha * k * beta(k - 1/c, 1 + 1/c);
else
    m = alpha * (2^(1/c) - 1)^(1/k);
end
end

% =====================================================================
%   ADAPTIVE BINS WITH t0 & PRE-T0 GUARANTEE
% =====================================================================
function [edges, bin] = make_bins_adaptive(t, target_count, max_bins, min_bins, t0_est, pre_t0_bins)
if nargin < 4 || isempty(min_bins),    min_bins    = 0;   end
if nargin < 5,                         t0_est      = [];  end
if nargin < 6 || isempty(pre_t0_bins), pre_t0_bins = 0;   end

t = t(:);
if isempty(t) || any(~isfinite(t)) || range(t)==0
    edges = [min(t); max(t)+eps];
    bin   = ones(size(t));
    return;
end

[ts,~] = sort(t);
N = numel(ts);

starts = 1:target_count:N;
starts = starts(starts <= N);

if numel(starts) > 2
    inner_idx   = starts(2:end-1);
    inner_edges = ts(inner_idx);
else
    inner_edges = [];
end

edges = unique([min(t); inner_edges(:); max(t)+eps], 'sorted');

while numel(edges)-1 > max_bins
    if numel(edges) <= 3, break; end
    keep = true(size(edges));
    interior = 2:numel(edges)-1;
    keep(interior(2:2:end)) = false;
    edges = edges(keep);
    edges(1)   = min(t);
    edges(end) = max(t)+eps;
end

if ~isempty(t0_est) && isfinite(t0_est)
    t0c  = min(max(t0_est, min(t)), max(t));
    edges = unique([edges(:); t0c], 'sorted');
end

while numel(edges)-1 < min_bins
    widths = diff(edges);
    [~,i] = max(widths);
    mid = (edges(i) + edges(i+1))/2;
    if mid > edges(i) && mid < edges(i+1)
        edges = sort([edges; mid]);
    else
        break;
    end
end

if ~isempty(t0_est) && isfinite(t0_est) && pre_t0_bins > 0
    t0c = min(max(t0_est, min(t)), max(t));
    right_edges = edges(2:end);
    pre_count   = sum(right_edges <= t0c);
    guard = 0;
    while pre_count < pre_t0_bins
        guard = guard + 1; if guard > 2000, break; end
        idx = find(right_edges <= t0c, 1, 'last');
        if isempty(idx) || idx < 1
            i = 1;
        else
            i = idx;
        end
        mid = (edges(i) + edges(i+1))/2;
        if ~(mid > edges(i) && mid < edges(i+1))
            mid = edges(i) + (edges(i+1)-edges(i))/3;
        end
        edges = sort([edges; mid]);
        right_edges = edges(2:end);
        pre_count   = sum(right_edges <= t0c);
        if numel(edges)-1 > max_bins + min_bins + 500, break; end
    end
end

edges = unique(edges(:).', 'sorted');
edges = edges([true diff(edges) > 0]);
edges(1)   = min(t);
edges(end) = max(t) + eps;

bin = discretize(t, edges);
end


% =====================================================================
%     ENFORCE MIN COUNT NEAR t0 & OPTIONAL SHRINK AT t0 EDGE
% =====================================================================
function edges = enforce_min_count_near_t0(edges, t, t0, min_count)
if isempty(t0) || ~isfinite(t0) || min_count<=0, return; end
[~, idx_t0] = min(abs(edges - t0));
if ~isfinite(edges(idx_t0)) || abs(edges(idx_t0) - t0) > 10*eps, return; end
if idx_t0 >= numel(edges), return; end

guard = 0;
while true
    guard = guard + 1; if guard > 1000, break; end
    bin_right = t > edges(idx_t0) & t <= edges(idx_t0+1);
    nright = sum(bin_right);
    if nright >= min_count, break; end

    if idx_t0+2 <= numel(edges)
        edges(idx_t0+1) = [];  % merge right
    elseif idx_t0-1 >= 1
        edges(idx_t0) = [];    % merge left, then reinsert t0
        edges = unique([edges(:); t0], 'sorted');
        [~, idx_t0] = min(abs(edges - t0));
    else
        break;
    end
end
end

function p_hat = shrink_p_near_t0(t_lo, t_hi, p_hat, n_i, t0, prior_strength)
if prior_strength <= 0 || isempty(t0) || ~isfinite(t0), return; end
touch = (abs(t_lo - t0) < 10*eps) | (abs(t_hi - t0) < 10*eps);
idx = find(touch);
for k = idx(:).'
    left  = k-1 >= 1;  right = k+1 <= numel(p_hat);
    if left && right
        p0 = mean([p_hat(k-1), p_hat(k+1)], 'omitnan');
    elseif left
        p0 = p_hat(k-1);
    elseif right
        p0 = p_hat(k+1);
    else
        p0 = mean(p_hat, 'omitnan');
    end
    if ~isfinite(p0), p0 = 0.5; end
    a0 = max(p0, 1e-6) * prior_strength;
    b0 = max(1-p0, 1e-6) * prior_strength;

    yk = max(min(p_hat(k)*n_i(k), n_i(k)), 0);  % pseudo-successes
    p_hat(k) = (a0 + yk) / (a0 + b0 + n_i(k));
end
end


% =====================================================================
%              IN-BIN EM FOR MIXTURE WEIGHT p (μ,σ fixed)
% =====================================================================
% function p = mix_weight_EM(x, mu1, mu2, sigma, p, iters)
% if nargin<6, iters=8; end
% n1 = normpdf(x,mu1,sigma);
% n2 = normpdf(x,mu2,sigma);
% 
% for it=1:iters
%     r = p.*n2 ./ ( (1-p).*n1 + p.*n2 + realmin );
%     p = mean(r);
%     p = min(max(p,1e-6),1-1e-6);
% end
% end
function p = mix_weight_EM(x, params1, params2, p, iters)
if nargin<5, iters=8; end

% Burr distribution PDFs
% params1 and params2 are each [alpha, c, k] from the Burr distribution
n1 = pdf('Burr', x, params1(1), params1(2), params1(3));
n2 = pdf('Burr', x, params2(1), params2(2), params2(3));

for it=1:iters
    r = p.*n2 ./ ( (1-p).*n1 + p.*n2 + realmin );
    p = mean(r);
    p = min(max(p,1e-6),1-1e-6);
end
end

% =====================================================================
%          CHOOSE TARGET COUNT PER SAMPLE (derived or fixed)
% =====================================================================
function tc = choose_target_count_per_sample(N, opts)
if strcmpi(opts.target_count_mode, 'by-total')
    tc = max( round(opts.target_count_scale * max(N,1)), opts.tc_min );
else
    tc = opts.target_count;
end
tc = max(1, tc);
end


% =====================================================================
%                    LOGISTIC 4P (weighted β-optimized)
% =====================================================================
function out = fit_bins_logistic4_optbeta(t, p, n, dispMode, wantCI, alpha, ...
                                          weight_mode, beta_single, ...
                                          do_beta_opt, beta_grid, weight_floor)
if nargin<5 || isempty(wantCI), wantCI = true; end
if nargin<6 || isempty(alpha),  alpha  = 0.05; end
if nargin<7, weight_mode = ''; end
if nargin<8 || isempty(beta_single), beta_single = 1; end
if nargin<9 || isempty(do_beta_opt), do_beta_opt = true; end
if nargin<10, beta_grid = []; end
if nargin<11 || isempty(weight_floor), weight_floor = 0; end

y = p .* n;

% init
B0   = max(min(p), 0);
T0   = min(max(p), 0.999);
A0   = max(T0 - B0, 1e-3);
t50_0= median(t);
k0   = 1 / max(iqr(t), eps);

theta0 = [B0; log(A0); t50_0; log(k0)];
nll0   = @(th) nll_logistic4_unw(th, t, y, n);

opts = optimset('Display', dispMode, 'MaxIter', 3000, 'MaxFunEvals', 6e4, ...
                'TolX',1e-8,'TolFun',1e-8);
[th0, ~] = fminsearch(nll0, theta0, opts);

B   = th0(1);  A=exp(th0(2));  t50=th0(3);  k=exp(th0(4));
[p_hat0, ll0, AIC0, BIC0, R20] = eval_fit_logistic(t, p, n, B, A, t50, k, th0);

final = struct('B',B,'A',A,'t50',t50,'k',k, ...
               'll',ll0,'AIC',AIC0,'BIC',BIC0,'R2',R20, ...
               'th',th0,'p_hat',p_hat0,'weights',ones(size(t)), ...
               'beta',0,'nll_fun',nll0);

% beta search ONLY over provided grid
if ~isempty(weight_mode) && ischar(weight_mode) && ~strcmpi(weight_mode,'') && do_beta_opt
    betas = beta_grid(:).';
    if isempty(betas), error('logistic_optimize_beta=true but logistic_beta_grid is empty.'); end
    best = final; best.R2 = -Inf; best.ll = -Inf;

    for bb = betas
        bb = max(bb,0);
        if bb==0
            cand = final; cand.beta = 0; cand.weights = ones(size(t));
        else
            w = compute_transition_weights_logistic(weight_mode, bb, p, t, B, A, t50, k, weight_floor);
            thetaW0 = [B; log(max(A,1e-12)); t50; log(max(k,1e-12))];
            nllW = @(th) nll_logistic4_w(th, t, y, n, w);
            [thW, ~] = fminsearch(nllW, thetaW0, opts);
            Bw  = thW(1); Aw=exp(thW(2)); t50w=thW(3); kw=exp(thW(4));
            [p_hatW, llW, AICW, BICW, R2W] = eval_fit_logistic(t, p, n, Bw, Aw, t50w, kw, thW);
            cand = struct('B',Bw,'A',Aw,'t50',t50w,'k',kw, ...
                          'll',llW,'AIC',AICW,'BIC',BICW,'R2',R2W, ...
                          'th',thW,'p_hat',p_hatW,'weights',w, ...
                          'beta',bb,'nll_fun',nllW);
        end

        if (cand.R2 > best.R2) || (abs(cand.R2 - best.R2) < 1e-12 && cand.ll > best.ll)
            best = cand;
        end
    end
    final = best;

elseif ~isempty(weight_mode) && ischar(weight_mode) && ~strcmpi(weight_mode,'') && ~do_beta_opt
    bb = max(beta_single,0);
    if bb==0
        % already final
    else
        w = compute_transition_weights_logistic(weight_mode, bb, p, t, B, A, t50, k, weight_floor);
        thetaW0 = [B; log(max(A,1e-12)); t50; log(max(k,1e-12))];
        nllW = @(th) nll_logistic4_w(th, t, y, n, w);
        [thW, ~] = fminsearch(nllW, thetaW0, opts);
        Bw  = thW(1); Aw=exp(thW(2)); t50w=thW(3); kw=exp(thW(4));
        [p_hatW, llW, AICW, BICW, R2W] = eval_fit_logistic(t, p, n, Bw, Aw, t50w, kw, thW);
        final = struct('B',Bw,'A',Aw,'t50',t50w,'k',kw, ...
                       'll',llW,'AIC',AICW,'BIC',BICW,'R2',R2W, ...
                       'th',thW,'p_hat',p_hatW,'weights',w, ...
                       'beta',bb,'nll_fun',nllW);
    end
end

out = struct();
out.model  = 'logistic4';
out.params = struct('B',final.B,'A',final.A,'t50',final.t50,'k',final.k, ...
                    'bottom',final.B,'top',final.B+final.A,'span',final.A);
out.loglik = final.ll; out.aic=final.AIC; out.bic=final.BIC; out.R2=final.R2;
out.fit    = struct('t',t,'p',p,'n',n,'p_hat',final.p_hat,'weights',final.weights);
out.notes  = struct('weight_mode',weight_mode, ...
                    'selected_beta',final.beta, ...
                    'beta_grid',ternary(do_beta_opt, beta_grid(:).', []), ...
                    'beta_optimized', do_beta_opt && ~isempty(weight_mode) && ~strcmpi(weight_mode,''));

% Uncertainty
if wantCI
    theta_hat = [final.B; log(max(final.A,1e-12)); final.t50; log(max(final.k,1e-12))];
    H  = numhess(final.nll_fun, theta_hat);
    H  = (H+H.')/2; [V,D]=eig(H); D=max(D,1e-9*eye(size(D))); H=V*D*V';
    Cov_theta = inv(H);
    A = final.A; K = final.k;
    J = diag([1, A, 1, K]);
    Cov_phi = J * Cov_theta * J.';
    SE = sqrt(max(diag(Cov_phi), 0));
    z = norminv(1-alpha/2);

    B = final.B; t50 = final.t50;
    ci = [ [B;A;t50;K] - z*SE, [B;A;t50;K] + z*SE ];
    out.Cov = Cov_phi;
    out.se  = struct('B',SE(1),'A',SE(2),'t50',SE(3),'k',SE(4), ...
                     'top',sqrt(SE(1)^2 + SE(2)^2), ...
                     'span',SE(2));
    out.ci  = struct('B',ci(1,:), 'A',ci(2,:), 't50',ci(3,:), 'k',ci(4,:), ...
                     'top',[ (B+A)-z*out.se.top, (B+A)+z*out.se.top ], ...
                     'span',[ A - z*out.se.span,  A + z*out.se.span ]);
end
end

function [p_hat, ll, AIC, BIC, R2] = eval_fit_logistic(t, p, n, B, A, t50, k, th)
p_hat = logistic4_eval(t, B, A, t50, k);
ll = sum( (p.*n) .* log(max(p_hat,eps)) + (n - p.*n) .* log(max(1-p_hat,eps)) );
pcount = numel(th);
AIC = 2*pcount - 2*ll; BIC = pcount*log(sum(n)) - 2*ll;
R2 = corr(p, p_hat, 'rows', 'complete')^2;
end

function v = nll_logistic4_unw(th, t, y, n)
B   = th(1); A=exp(th(2)); t50=th(3); k=exp(th(4));
p   = logistic4_eval(t, B, A, t50, k);
v   = -sum( y .* log(max(p,eps)) + (n-y) .* log(max(1-p,eps)) );
end

function v = nll_logistic4_w(th, t, y, n, w)
B   = th(1); A=exp(th(2)); t50=th(3); k=exp(th(4));
p   = logistic4_eval(t, B, A, t50, k);
v   = -sum( w .* ( y .* log(max(p,eps)) + (n-y) .* log(max(1-p,eps)) ) );
end

function w = compute_transition_weights_logistic(mode, beta, p_obs, t, B, A, t50, k, wfloor)
if nargin<9 || isempty(wfloor), wfloor = 0; end
beta = max(beta,0);
Aeff = max(A, 1e-9);

switch lower(string(mode))
    case "p-mid"
        q = (p_obs - B) ./ Aeff;
        q = min(max(q, 0), 1);
        w = (q .* (1 - q)) .^ beta;

    case "p-mid-model"
        s  = 1 ./ (1 + exp(-k .* (t - t50)));
        p_hat = B + A .* s;
        q = (p_hat - B) ./ Aeff;
        q = min(max(q, 0), 1);
        w = (q .* (1 - q)) .^ beta;

    case "slope"
        s  = 1 ./ (1 + exp(-k .* (t - t50)));
        dp = A .* k .* s .* (1 - s);
        w  = abs(dp) .^ beta;

    otherwise
        w = ones(size(t));
end

w(~isfinite(w)) = 0;
m = mean(w(w>0));
if ~isfinite(m) || m<=0, m = 1; end
w = w ./ m;

if wfloor > 0
    w = max(w, wfloor);
    m2 = mean(w);
    if isfinite(m2) && m2>0, w = w ./ m2; end
end
end


% =====================================================================
%     GAMMA INIT FROM LOGISTIC VIA QUANTILE MATCHING (inner fit)
% =====================================================================
function init = init_gamma_from_logistic_quintiles(t, B, A, t50, k, t0_fixed, quants, wq, dispMode)
% Build target times for p=B+q*A from logistic invert, then fit (N,k)
% (t0 fixed if provided) to match gamcdf(t - t0; N, 1/k) ≈ q at those times.

if nargin<7 || isempty(quants), quants = [0.2 0.4 0.6 0.8]; end
quants = quants(:).';
if nargin<8 || isempty(wq), wq = ones(size(quants)); end
wq = wq(:).'/sum(wq);

% Inverse logistic for times at each q
tq = zeros(size(quants));
for i=1:numel(quants)
    q = min(max(quants(i),1e-6), 1-1e-6);
    tq(i) = logistic_inverse_time(B, A, t50, k, q);
end

% Centering on t0
if ~isempty(t0_fixed) && isfinite(t0_fixed)
    tau = max(tq - t0_fixed, 0);
else
    % If t0 is unknown here, estimate a rough t0 as the earliest tq - small epsilon
    t0_guess = min(tq) - 0.25*iqr(t); % rough
    tau = max(tq - t0_guess, 0);
end

% Initial guess for N,k using moments of logistic slope width:
wL = iqr(t) / max(k,1e-9);    % crude width proxy
N0 = 2 + 2*(range(t)/max(wL,eps));     % shape grows with narrower transition
N0 = max(1.2, min(40, N0));
k0 = 1 / max(wL, eps);                    % rate increases with narrower width

theta0 = [log(N0); log(k0)];

if ~isempty(t0_fixed) && isfinite(t0_fixed)
    nll = @(th) sum( wq .* ( gamcdf(tq - t0_fixed,  max(exp(th(1)),1e-6), 1/max(exp(th(2)),1e-9)) - quants ).^2 );
else
    % jointly fit t0 too (mildly) for the inner init, then pass results forward
    theta0 = [theta0; median(tq)-0.5*iqr(t)];
    nll = @(th) sum( wq .* ( gamcdf(tq - th(3), max(exp(th(1)),1e-6), 1/max(exp(th(2)),1e-9)) - quants ).^2 );
end

opt = optimset('Display', dispMode, 'MaxIter', 2000, 'MaxFunEvals', 20000, ...
               'TolX',1e-8, 'TolFun',1e-8);
[th, ~] = fminsearch(nll, theta0, opt);

if ~isempty(t0_fixed) && isfinite(t0_fixed)
    N_init = max(exp(th(1)),1e-6);
    k_init = max(exp(th(2)),1e-9);
    t0_init = t0_fixed;
else
    N_init = max(exp(th(1)),1e-6);
    k_init = max(exp(th(2)),1e-9);
    t0_init = th(3);
end

init = struct('B',B,'A',A,'t0',t0_init,'N',N_init,'k',k_init);
end

function t_q = logistic_inverse_time(B, A, t50, k, q)
q = min(max(q, 1e-9), 1-1e-9);
% q here is fraction of span: p = B + q*A
% Logistic: p = B + A / (1 + exp(-k (t - t50))) -> solve for t
r = A ./ (q*A);  % = 1/q
% 1 + exp(-kΔ) = r -> exp(-kΔ)= r-1 -> Δ = -log(r-1)/k
t_q = t50 - (1./max(k,1e-12)) .* log(max(r-1, 1e-12));
end


% =====================================================================
%            GAMMA + ONSET (weighted β-optimized, t0 optional fixed)
% =====================================================================
function out = fit_bins_gamma_onset_optbeta(t, p, n, dispMode, wantCI, ...
                                            t0_info, t0_se_info, ...
                                            t0_mode, prior_scale, se_min, se_max, when_missing, ...
                                            weight_mode, beta_single, ...
                                            do_beta_opt, beta_grid, weight_floor)
% Prior-aware gamma fit with optional slope/p-mid weighting and per-sample beta search.
%
% t0_mode: 'fixed' | 'free' | 'prior'
%   'fixed' -> t0 locked to t0_info
%   'free'  -> t0 estimated with no constraint
%   'prior' -> t0 estimated with Gaussian prior N(t0_info, (se_eff)^2), where se_eff is
%              clamped to [se_min,se_max] and the penalty is scaled by prior_scale.
%
% Notes:
%  * The prior only affects the objective (negative log-posterior) term; data weights (w)
%    still apply only to the data likelihood contribution.
%  * If t0_se_info is missing/NaN and t0_mode='prior', we switch to when_missing ('free' or 'fixed').

if nargin<5 || isempty(wantCI), wantCI = true; end
if nargin<7, t0_se_info = []; end
if nargin<8 || isempty(t0_mode), t0_mode = 'fixed'; end
if nargin<9 || isempty(prior_scale), prior_scale = 1.0; end
if nargin<10 || isempty(se_min), se_min = 0.05; end
if nargin<11 || isempty(se_max), se_max = Inf; end
if nargin<12 || isempty(when_missing), when_missing = 'free'; end
if nargin<13, weight_mode = ''; end
if nargin<14 || isempty(beta_single), beta_single = 1; end
if nargin<15 || isempty(do_beta_opt), do_beta_opt = true; end
if nargin<16, beta_grid = []; end
if nargin<17 || isempty(weight_floor), weight_floor = 0; end

y = p .* n;

% --- inits
B0  = max(min(p), 0);
T0  = min(max(p), 0.999);
A0  = max(T0 - B0, 1e-3);
[~,ord] = sort(t); ts=t(ord); ps=p(ord);
thr  = B0 + 0.1*A0; idx = find(ps > thr, 1, 'first');
t0_0 = ternary(~isempty(idx), ts(max(1,idx-1)), min(t));
N0   = 2;
k0   = 1 / max(4*iqr(t), eps);

% --- decide mode & effective prior
t0_has = ~isempty(t0_info) && isfinite(t0_info);
se_has = ~isempty(t0_se_info) && isfinite(t0_se_info) && (t0_se_info>0);
use_fixed = strcmpi(t0_mode,'fixed') && t0_has;
use_prior = strcmpi(t0_mode,'prior') && t0_has && se_has;
if strcmpi(t0_mode,'prior') && (~t0_has || ~se_has)
    % fall back if we can't form a prior
    if strcmpi(when_missing,'fixed') && t0_has
        use_fixed = true; use_prior = false;
    else
        use_fixed = false; use_prior = false; % -> free
    end
end

% clamp SE and scale
if use_prior
    se_eff = min(max(t0_se_info, se_min), se_max);
    lam = prior_scale;  % penalty multiplier (>=0)
else
    se_eff = NaN; lam = 0;
end

% --- baseline objective & start
if use_fixed
    theta0 = [B0; log(A0); log(N0); log(k0)];
    nll0 = @(th) nll_gamma_onset_fixed(th, t, y, n, t0_info);
elseif use_prior
    theta0 = [B0; log(A0); t0_info; log(N0); log(k0)];
    nll0 = @(th) nll_gamma_onset_prior(th, t, y, n, t0_info, se_eff, lam);
else % free
    theta0 = [B0; log(A0); t0_0;   log(N0); log(k0)];
    nll0 = @(th) nll_gamma_onset(th, t, y, n);
end

opts = optimset('Display', dispMode, 'MaxIter', 4000, 'MaxFunEvals', 8e4, ...
                'TolX',1e-8,'TolFun',1e-8);
[th0, ~] = fminsearch(nll0, theta0, opts);

if use_fixed
    B  = th0(1); A=exp(th0(2)); t0=t0_info; N=max(exp(th0(3)),1e-6); k=exp(th0(4));
else
    B  = th0(1); A=exp(th0(2)); t0=th0(3);    N=max(exp(th0(4)),1e-6); k=exp(th0(5));
end

[p_hat0, ll0, AIC0, BIC0, R20] = eval_fit_gamma(t, p, n, B, A, t0, N, k, th0);

final = struct('B',NaN,'A',NaN,'t0',NaN,'N',NaN,'k',NaN, ...
               'll',-Inf,'AIC',Inf,'BIC',Inf,'R2',-Inf, ...
               'th',[],'p_hat',[],'weights',[],'beta',[], ...
               'nll_fun',[]);

% --- beta search ONLY over provided grid if weighting is requested
if ~isempty(weight_mode) && ischar(weight_mode) && ~strcmpi(weight_mode,'') && do_beta_opt
    betas = beta_grid(:).';
    if isempty(betas)
        error('gamma_optimize_beta=true but gamma_beta_grid is empty.');
    end

    for bb = betas
        bb = max(bb,0);
        if bb==0
            cand = struct('B',B,'A',A,'t0',t0,'N',N,'k',k, ...
                          'll',ll0,'AIC',AIC0,'BIC',BIC0,'R2',R20, ...
                          'th',th0,'p_hat',p_hat0,'weights',ones(size(t)), ...
                          'beta',0,'nll_fun',nll0);
        else
            % (we keep your existing weighting function — slope/p-mid/etc.)
            w = compute_transition_weights(weight_mode, bb, p, t, B, A, t0, N, k, weight_floor);

            if use_fixed
                thetaW0 = [B; log(max(A,1e-12)); log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_fixed_w(th, t, y, n, t0_info, w);
            elseif use_prior
                thetaW0 = [B; log(max(A,1e-12)); t0; log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_prior_w(th, t, y, n, w, t0_info, se_eff, lam);
            else
                thetaW0 = [B; log(max(A,1e-12)); t0; log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_w(th, t, y, n, w);
            end

            [thW, ~] = fminsearch(nllW, thetaW0, opts);
            if use_fixed
                Bw  = thW(1); Aw=exp(thW(2)); t0w=t0_info; Nw=max(exp(thW(3)),1e-6); kw=exp(thW(4));
            else
                Bw  = thW(1); Aw=exp(thW(2)); t0w=thW(3);   Nw=max(exp(thW(4)),1e-6); kw=exp(thW(5));
            end

            [p_hatW, llW, AICW, BICW, R2W] = eval_fit_gamma(t, p, n, Bw, Aw, t0w, Nw, kw, thW);
            cand = struct('B',Bw,'A',Aw,'t0',t0w,'N',Nw,'k',kw, ...
                          'll',llW,'AIC',AICW,'BIC',BICW,'R2',R2W, ...
                          'th',thW,'p_hat',p_hatW,'weights',w,'beta',bb, ...
                          'nll_fun',nllW);
        end

        if (cand.R2 > final.R2) || (abs(cand.R2 - final.R2) < 1e-12 && cand.ll > final.ll)
            final = cand;
        end
    end

else
    % no optimization (single beta or unweighted)
    if isempty(weight_mode) || strcmpi(weight_mode,'')
        final = struct('B',B,'A',A,'t0',t0,'N',N,'k',k, ...
                       'll',ll0,'AIC',AIC0,'BIC',BIC0,'R2',R20, ...
                       'th',th0,'p_hat',p_hat0,'weights',[], ...
                       'beta',0,'nll_fun',nll0);
    else
        bb = max(beta_single,0);
        if bb==0
            final = struct('B',B,'A',A,'t0',t0,'N',N,'k',k, ...
                           'll',ll0,'AIC',AIC0,'BIC',BIC0,'R2',R20, ...
                           'th',th0,'p_hat',p_hat0,'weights',ones(size(t)), ...
                           'beta',0,'nll_fun',nll0);
        else
            w = compute_transition_weights(weight_mode, bb, p, t, B, A, t0, N, k, weight_floor);
            if use_fixed
                thetaW0 = [B; log(max(A,1e-12)); log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_fixed_w(th, t, y, n, t0_info, w);
            elseif use_prior
                thetaW0 = [B; log(max(A,1e-12)); t0; log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_prior_w(th, t, y, n, w, t0_info, se_eff, lam);
            else
                thetaW0 = [B; log(max(A,1e-12)); t0; log(max(N,1e-12)); log(max(k,1e-12))];
                nllW = @(th) nll_gamma_onset_w(th, t, y, n, w);
            end
            [thW, ~] = fminsearch(nllW, thetaW0, opts);
            if use_fixed
                Bw  = thW(1); Aw=exp(thW(2)); t0w=t0_info; Nw=max(exp(thW(3)),1e-6); kw=exp(thW(4));
            else
                Bw  = thW(1); Aw=exp(thW(2)); t0w=thW(3);   Nw=max(exp(thW(4)),1e-6); kw=exp(thW(5));
            end
            [p_hatW, llW, AICW, BICW, R2W] = eval_fit_gamma(t, p, n, Bw, Aw, t0w, Nw, kw, thW);
            final = struct('B',Bw,'A',Aw,'t0',t0w,'N',Nw,'k',kw, ...
                           'll',llW,'AIC',AICW,'BIC',BICW,'R2',R2W, ...
                           'th',thW,'p_hat',p_hatW,'weights',w,'beta',bb, ...
                           'nll_fun',nllW);
        end
    end
end

% --- pack outputs
out = struct();
out.model  = 'gamma';
out.params = struct('B',final.B,'A',final.A,'t0',final.t0,'N',final.N,'k',final.k, ...
                    'bottom',final.B,'top',final.B+final.A,'span',final.A);
out.loglik = final.ll; out.aic=final.AIC; out.bic=final.BIC; out.R2=final.R2;
out.fit    = struct('t',t,'p',p,'n',n,'p_hat',final.p_hat,'weights',final.weights);
out.notes  = struct('t0_mode', t0_mode, ...
                    't0_used', ternary(use_fixed, t0_info, final.t0), ...
                    't0_prior_used', use_prior, ...
                    't0_prior_center', ternary(use_prior, t0_info, NaN), ...
                    't0_prior_se_eff', ternary(use_prior, se_eff, NaN), ...
                    't0_prior_scale',   ternary(use_prior, prior_scale, NaN), ...
                    'weight_mode',weight_mode, ...
                    'selected_beta',final.beta, ...
                    'beta_grid',ternary(do_beta_opt, beta_grid(:).', []), ...
                    'beta_optimized', do_beta_opt && ~isempty(weight_mode) && ~strcmpi(weight_mode,''));

% --- uncertainty (Hessian includes prior if used)
if wantCI
    z = 1.96;
    if use_fixed
        theta_hat = [final.B; log(max(final.A,1e-12)); log(max(final.N,1e-12)); log(max(final.k,1e-12))];
    else
        theta_hat = [final.B; log(max(final.A,1e-12)); final.t0; log(max(final.N,1e-12)); log(max(final.k,1e-12))];
    end

    H  = numhess(final.nll_fun, theta_hat);
    H  = (H+H.')/2; [V,D]=eig(H); D=max(D,1e-9*eye(size(D))); H=V*D*V';
    Cov_theta = inv(H);

    if use_fixed
        A=final.A; N=final.N; k=final.k;
        J = diag([1, A, N, k]);                 % [B logA logN logk] -> [B A N k]
        Cov_phi = J * Cov_theta * J.'; SE = sqrt(max(diag(Cov_phi),0));
        out.Cov = Cov_phi;
        out.se  = struct('B',SE(1),'A',SE(2),'t0',0,'N',SE(3),'k',SE(4), ...
                         'top',sqrt(SE(1)^2 + SE(2)^2),'span',SE(2));
        B=final.B; A=final.A; t0=final.t0;
        out.ci  = struct('B',[B - z*SE(1), B + z*SE(1)], ...
                         'A',[A - z*SE(2), A + z*SE(2)], ...
                         't0',[t0 t0], ...
                         'N',[final.N - z*SE(3), final.N + z*SE(3)], ...
                         'k',[final.k - z*SE(4), final.k + z*SE(4)], ...
                         'top',[ (B+A) - z*out.se.top, (B+A) + z*out.se.top ], ...
                         'span',[ A - z*out.se.span,    A + z*out.se.span ]);
    else
        A=final.A; N=final.N; k=final.k;
        J = diag([1, A, 1, N, k]);              % [B logA t0 logN logk] -> [B A t0 N k]
        Cov_phi = J * Cov_theta * J.'; SE = sqrt(max(diag(Cov_phi),0));
        out.Cov = Cov_phi;
        out.se  = struct('B',SE(1),'A',SE(2),'t0',SE(3),'N',SE(4),'k',SE(5), ...
                         'top',sqrt(SE(1)^2 + SE(2)^2),'span',SE(2));
        B=final.B; A=final.A; t0=final.t0;
        out.ci  = struct('B',[B - z*SE(1), B + z*SE(1)], ...
                         'A',[A - z*SE(2), A + z*SE(2)], ...
                         't0',[t0 - z*SE(3), t0 + z*SE(3)], ...
                         'N',[final.N - z*SE(4), final.N + z*SE(4)], ...
                         'k',[final.k - z*SE(5), final.k + z*SE(5)], ...
                         'top',[ (B+A) - z*out.se.top, (B+A) + z*out.se.top ], ...
                         'span',[ A - z*out.se.span,    A + z*out.se.span ]);
    end
end
end

function [p_hat, ll, AIC, BIC, R2] = eval_fit_gamma(t, p, n, B, A, t0, N, k, th)
% EVAL_FIT_GAMMA
% Evaluate likelihood & simple GOF metrics for the onset-gamma CDF model.
% Inputs:
%   t,p,n : binned time, fraction, counts
%   B,A,t0,N,k : gamma-onset parameters
%   th     : parameter vector used during the fit (for counting params)
% Outputs:
%   p_hat : model-predicted p(t)
%   ll    : binomial log-likelihood
%   AIC,BIC: standard information criteria
%   R2    : squared Pearson correlation between p and p_hat

p_hat = gamma_onset_eval(t, B, A, t0, N, k);

% Binomial log-likelihood
y = p .* n;
ll = sum( y .* log(max(p_hat,eps)) + (n - y) .* log(max(1 - p_hat, eps)) );

% ICs use number of *free* params in 'th'
pcount = numel(th);
AIC = 2*pcount - 2*ll;
BIC = pcount*log(sum(n)) - 2*ll;

% Simple R^2 on proportions
R2 = corr(p, p_hat, 'rows', 'complete')^2;
end


% ---- NLLs with/without prior and with/without weights ----

function v = nll_gamma_onset(th, t, y, n)
B=th(1); A=exp(th(2)); t0=th(3);
N=max(exp(th(4)),1e-6); k=exp(th(5));
p=gamma_onset_eval(t,B,A,t0,N,k);
v=-sum( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) );
end

function v = nll_gamma_onset_w(th, t, y, n, w)
B=th(1); A=exp(th(2)); t0=th(3);
N=max(exp(th(4)),1e-6); k=exp(th(5));
p=gamma_onset_eval(t,B,A,t0,N,k);
v=-sum( w .* ( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) ) );
end

function v = nll_gamma_onset_fixed(th, t, y, n, t0_fixed)
B=th(1); A=exp(th(2));
N=max(exp(th(3)),1e-6); k=exp(th(4));
p=gamma_onset_eval(t,B,A,t0_fixed,N,k);
v=-sum( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) );
end

function v = nll_gamma_onset_fixed_w(th, t, y, n, t0_fixed, w)
B=th(1); A=exp(th(2));
N=max(exp(th(3)),1e-6); k=exp(th(4));
p=gamma_onset_eval(t,B,A,t0_fixed,N,k);
v=-sum( w .* ( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) ) );
end

% ---- Prior on t0: Gaussian N(t0c, se^2); lam scales the penalty strength
function v = nll_gamma_onset_prior(th, t, y, n, t0c, se, lam)
B=th(1); A=exp(th(2)); t0=th(3);
N=max(exp(th(4)),1e-6); k=exp(th(5));
p=gamma_onset_eval(t,B,A,t0,N,k);
nll = -sum( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) );
pen = 0.5 * lam * ((t0 - t0c)/max(se,eps))^2;
v = nll + pen;
end

function v = nll_gamma_onset_prior_w(th, t, y, n, w, t0c, se, lam)
B=th(1); A=exp(th(2)); t0=th(3);
N=max(exp(th(4)),1e-6); k=exp(th(5));
p=gamma_onset_eval(t,B,A,t0,N,k);
nll = -sum( w .* ( y.*log(max(p,eps)) + (n-y).*log(max(1-p,eps)) ) );
pen = 0.5 * lam * ((t0 - t0c)/max(se,eps))^2;
v = nll + pen;
end


% =====================================================================
%   GAMMA WEIGHTS (midness/slope in gamma space; 'slope' usually from L)
% =====================================================================
function w = compute_transition_weights_gamma(mode, beta, p_obs, t, B, A, t0, N, k, wfloor)
if nargin<10 || isempty(wfloor), wfloor = 0; end
beta = max(beta,0);
Aeff = max(A, 1e-9);

switch lower(string(mode))
    case "p-mid"
        q = (p_obs - B) ./ Aeff;
        q = min(max(q, 0), 1);
        w = (q .* (1 - q)) .^ beta;

    case "p-mid-model"
        tau = max(t - t0, 0);
        p_hat = B + A .* gamcdf(tau, max(N,1e-6), 1/max(k,1e-9));
        q = (p_hat - B) ./ Aeff;
        q = min(max(q, 0), 1);
        w = (q .* (1 - q)) .^ beta;

    case "slope"
        % If called, prefer using logistic-slope weights via logistic helper.
        % This function remains for completeness if 'slope' is asked but no L is available.
        tau = max(t - t0, 0);
        w = (gampdf(tau, max(N,1e-6), 1/max(k,1e-9))) .^ beta;

    otherwise
        w = ones(size(t));
end

w(~isfinite(w)) = 0;
m = mean(w(w>0));
if ~isfinite(m) || m<=0, m = 1; end
w = w ./ m;

if wfloor > 0
    w = max(w, wfloor);
    m2 = mean(w);
    if isfinite(m2) && m2>0, w = w ./ m2; end
end
end


% =====================================================================
%                       CURVE EVALUATORS
% =====================================================================
function p = logistic4_eval(t, B, A, t50, k)
p = B + A ./ (1 + exp(-k .* (t - t50)));
p = min(max(p, 1e-9), 1 - 1e-9);
end

function p = gamma_onset_eval(t,B,A,t0,N,k)
tau=max(t - t0,0);
theta=1 / max(k, 1e-9);
p=B + A .* gamcdf(tau, max(N,1e-6), theta);
p=min(max(p, 1e-9), 1-1e-9);
end


% ======================================================================
%             t0 detection from event rate (+ uncertainty)
% ======================================================================
function [t0_hat, se_t0] = estimate_t0_from_rate(t, binw, smooth_w, direction, search_window, do_boot, B, minNboot, alpha)
t = t(:);
t0_hat = median(t); se_t0 = NaN;
if isempty(t) || ~isfinite(min(t)) || range(t)==0
    return;
end
tmin = min(t); tmax = max(t);
if nargin>=5 && ~isempty(search_window)
    tmin = max(tmin, search_window(1));
    tmax = min(tmax, search_window(2));
end
if nargin<6 || isempty(do_boot), do_boot = true; end
if nargin<7 || isempty(B),       B = 200; end
if nargin<8 || isempty(minNboot),minNboot = 200; end
if nargin<9 || isempty(alpha),   alpha = 0.05; end

edges = tmin:binw:tmax; if numel(edges) < 5, edges = linspace(tmin,tmax,5); end
[counts,~,~] = histcounts(t, edges);
rate = counts(:) / binw;
rate_s = movmean(rate, max(1,smooth_w), 'Endpoints','shrink');
dr = [0; diff(rate_s)];
if strcmpi(direction,'decrease'), [~, idx] = min(dr); else, [~, idx] = max(dr); end
centers = (edges(1:end-1) + edges(2:end))/2;
idx = max(1, min(idx, numel(centers)));
t0_hat = centers(idx);

if do_boot && numel(t) >= minNboot
    t0s = nan(B,1);
    n = numel(t);
    for b=1:B
        tb = t(randi(n, n, 1));
        tbmin = min(tb); tbmax = max(tb);
        edges_b = tbmin:binw:tbmax; if numel(edges_b) < 5, edges_b = linspace(tbmin,tbmax,5); end
        [counts_b,~,~] = histcounts(tb, edges_b);
        rate_b = counts_b(:) / binw;
        rate_bs = movmean(rate_b, max(1,smooth_w), 'Endpoints','shrink');
        dr_b = [0; diff(rate_bs)];
        if strcmpi(direction,'decrease'), [~, idxb] = min(dr_b); else, [~, idxb] = max(dr_b); end
        centers_b = (edges_b(1:end-1) + edges_b(2:end))/2;
        if isempty(centers_b), t0s(b) = NaN; else
            idxb = max(1, min(idxb, numel(centers_b)));
            t0s(b) = centers_b(idxb);
        end
    end
    t0s = t0s(isfinite(t0s));
    if numel(t0s) >= max(10, 0.1*B)
        se_t0 = std(t0s, 0);
    else
        se_t0 = binw / sqrt(12);
    end
else
    se_t0 = binw / sqrt(12);
end
end


% =====================================================================
%                         UTILITIES
% =====================================================================
function H = numhess(fun, theta, epsv)
if nargin<3, epsv = 1e-5; end
theta = theta(:); n = numel(theta);
H = zeros(n); f0 = fun(theta);
h = epsv * max(1, abs(theta));
for i=1:n
    ei = zeros(n,1); ei(i)=h(i);
    fpi = fun(theta + ei); fmi = fun(theta - ei);
    H(i,i) = (fpi - 2*f0 + fmi)/(h(i)^2);
end
for i=1:n
    ei = zeros(n,1); ei(i)=h(i);
    for j=i+1:n
        ej = zeros(n,1); ej(j)=h(j);
        fpp = fun(theta + ei + ej);
        fpm = fun(theta + ei - ej);
        fmp = fun(theta - ei + ej);
        fmm = fun(theta - ei - ej);
        Hij = (fpp - fpm - fmp + fmm) / (4*h(i)*h(j));
        H(i,j)=Hij; H(j,i)=Hij;
    end
end
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
