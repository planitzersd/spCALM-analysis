function res = tq_isotonic_struct_2(S, opts)
% tq_isotonic_struct  Robust time-to-fraction (e.g., t50) from noisy, monotone-ish data
% with optional per-sample bootstrap CIs/SEs, using TIME-BASED baseline/plateau windows.
%
%   res = tq_isotonic_struct(S)
%   res = tq_isotonic_struct(S, opts)
%
% INPUT
%   S    : structure array (N samples) with fields per sample:
%          - t    : [ni x 1] bin centers (optional if t_lo/t_hi provided)
%          - r    : [ni x 1] response per bin
%          - n    : [ni x 1] (optional) counts per bin (weights)
%          - t_lo : [ni x 1] (optional) bin lower edges
%          - t_hi : [ni x 1] (optional) bin upper edges
%   opts : (optional) struct
%          % threshold
%          - q            : scalar or vector in (0,1); default 0.5 (t50)
%          - kSigma       : detectability multiplier for sigma (MAD); default 3
%          - useCountsW   : use S(i).n as weights in isotonic? default true
%          - minBins      : min bins to attempt fit; default 3
%          % TIME WINDOWS (seconds)
%          - base_sec     : baseline window length from START (seconds)  [default: []]
%          - plat_sec     : plateau  window length to END   (seconds)    [default: []]
%            If either is empty, falls back to fractions below.
%          % (fallback) FRACTIONS of bins (kept for backward compat.)
%          - baseFrac     : fraction earliest bins if base_sec empty; default 0.15
%          - platFrac     : fraction latest   bins if plat_sec empty; default 0.15
%          % bootstrap (per-sample)
%          - boot.B       : number of bootstrap reps; default 0 (off)
%          - boot.alpha   : CI level; default 0.05 (→ 95% CI)
%          - boot.seed    : rng seed (empty = leave RNG as-is)
%          - boot.scheme  : 'wild' (default) or 'pairs'
%
% OUTPUT
%   If opts.q is scalar: res is [N x 1] struct with fields
%     tq, censored, se, ci(1x2), cross_rate, b, p, A, sigma, theta, idxCross, nBins, q
%   If opts.q is vector: res is a struct with fields per q (e.g., res.q50),
%     each an [N x 1] struct array as above.
% -------------------------------------------------------------------------

    if nargin < 2, opts = struct(); end
    opts = filldefaults(opts, struct( ...
        'q',          0.5, ...
        'kSigma',     3.0, ...
        'useCountsW', true, ...
        'minBins',    3, ...
        'base_sec',   [], ...
        'plat_sec',   [], ...
        'baseFrac',   0.15, ...
        'platFrac',   0.15, ...
        'boot',       struct('B',500,'alpha',0.05,'seed',123, 'scheme','wild')));

    if ~isfield(opts,'boot') || isempty(opts.boot), opts.boot = struct(); end
    opts.boot = filldefaults(opts.boot, struct('B',0,'alpha',0.05,'seed',[], 'scheme','wild'));

    qvec = opts.q(:)';  % allow multiple thresholds
    N = numel(S);

    % Prototype result
    proto = struct('name', NaN, 'tq', NaN, 'censored', true, ...
                   'se', NaN, 'ci', [NaN NaN], 'cross_rate', NaN, ...
                   'b', NaN, 'p', NaN, 'A', NaN, 'sigma', NaN, 'theta', NaN, ...
                   'idxCross', NaN, 'nBins', NaN, 'q', NaN);

    out = struct();
    for q = qvec
        % Pre-allocate for valid samples only
        validIdx = false(N, 1);
        Rq_temp = repmat(proto, N, 1);
        
        if ~isempty(opts.boot.seed), rng(opts.boot.seed); end

        for i = 1:N
            % ---- Try to extract this sample (skip if invalid) ----
            try
                [t, r, w, tlo, thi] = extract_trw_edges(S(i), opts.useCountsW);
            catch ME
                % Sample has empty/invalid fields - skip it
                warning('tq_isotonic_struct:InvalidSample', ...
                    'Skipping sample %d: %s', i, ME.message);
                continue;
            end
            
            good = isfinite(t) & isfinite(r) & isfinite(w) & (w >= 0);
            t=t(good); r=r(good); w=w(good);
            if ~isempty(tlo), tlo=tlo(good); end
            if ~isempty(thi), thi=thi(good); end
            Rq_temp(i).nBins = numel(t); Rq_temp(i).q = q;

            if numel(t) < opts.minBins, continue; end
            
            % Mark this sample as valid
            validIdx(i) = true;

            % sort by time (and reorder edges)
            [t, ord] = sort(t(:));
            r = r(ord); w = w(ord);
            if ~isempty(tlo), tlo = tlo(ord); end
            if ~isempty(thi), thi = thi(ord); end
            if all(w==0), w(:) = 1; else, w(w==0) = eps; end

            % compute per-bin widths (sec)
            widths = bin_widths(t, tlo, thi);

            % ---- choose baseline/plateau bins by TIME windows ----
            [idxB, idxP] = select_time_windows(t, widths, opts.base_sec, opts.plat_sec, opts.baseFrac, opts.platFrac);

            % ---- isotonic fit ----
            yhat = weighted_pav(r, w);

            % ---- baseline/plateau & noise ----
            b  = wmedian(yhat(idxB), w(idxB));
            p  = wmedian(yhat(idxP), w(idxP));
            A  = p - b;

            sigma = mad(r(idxB), 0);  % robust baseline noise from RAW values
            if ~isfinite(sigma) || sigma==0, sigma = eps; end

            theta = b + q*A;

            % ---- point estimate (and censoring) ----
            [tq_hat, cens_hat, idxCross] = crossing_time_from_fit(t, yhat, theta);
            if A < opts.kSigma * sigma
                tq_hat = t(end); cens_hat = true; idxCross = NaN;
            end
            
            if isfield(S(i), 'name')
                Rq_temp(i).name = S(i).name;
            end
            Rq_temp(i).tq       = tq_hat;
            Rq_temp(i).censored = cens_hat;
            Rq_temp(i).b        = b;
            Rq_temp(i).p        = p;
            Rq_temp(i).A        = A;
            Rq_temp(i).sigma    = sigma;
            Rq_temp(i).theta    = theta;
            Rq_temp(i).idxCross = idxCross;

            % ---- bootstrap uncertainty (optional) ----
            B = opts.boot.B;
            if B > 0
                switch lower(opts.boot.scheme)
                    case 'wild'
                        [se, ci, cr] = bootstrap_wild_timewin(t, r, w, yhat, q, sigma, opts.kSigma, B, opts.boot.alpha, ...
                                                              tlo, thi, opts.base_sec, opts.plat_sec, opts.baseFrac, opts.platFrac);
                    case 'pairs'
                        [se, ci, cr] = bootstrap_pairs_timewin(t, r, w, q, sigma, opts.kSigma, B, opts.boot.alpha, ...
                                                               tlo, thi, opts.base_sec, opts.plat_sec, opts.baseFrac, opts.platFrac);
                    otherwise
                        error('Unknown boot.scheme: %s (use ''wild'' or ''pairs'')', opts.boot.scheme);
                end
                Rq_temp(i).se = se;
                Rq_temp(i).ci = ci;
                Rq_temp(i).cross_rate = cr;
            end
        end

        % Keep only valid samples in output
        Rq = Rq_temp(validIdx);
        out.(qfieldname(q)) = Rq;
    end

    if numel(qvec) == 1
        res = out.(qfieldname(qvec));
    else
        res = out;
    end
end

% ----------------- Bootstrap helpers (time-window aware) -----------------

function [se, ci, cross_rate] = bootstrap_wild_timewin(t, r, w, yhat, q, sigma_base, kSigma, B, alpha, ...
                                                       tlo, thi, base_sec, plat_sec, baseFrac, platFrac)
    n = numel(t);
    res = r - yhat; t_last = t(end);
    tq = nan(B,1); crossed = false(B,1);

    for b = 1:B
        sgn = 2*(rand(n,1) < 0.5) - 1;       % Rademacher ±1
        r_b = yhat + res .* sgn;
        yhat_b = weighted_pav(r_b, w);

        widths_b = bin_widths(t, tlo, thi);
        [idxB, idxP] = select_time_windows(t, widths_b, base_sec, plat_sec, baseFrac, platFrac);

        b_b = wmedian(yhat_b(idxB), w(idxB));
        p_b = wmedian(yhat_b(idxP), w(idxP));
        A_b = p_b - b_b;
        theta_b = b_b + q*A_b;

        if A_b < kSigma * sigma_base
            tq(b) = t_last; crossed(b) = false;
        else
            [tb, cens] = crossing_time_from_fit(t, yhat_b, theta_b);
            tq(b) = tb; crossed(b) = ~cens;
        end
    end

    cross_rate = mean(crossed);
    if any(crossed)
        x = tq(crossed);
        se = std(x);
        lo = quantile(x, alpha/2);
        hi = quantile(x, 1 - alpha/2);
        ci = [lo, hi];
    else
        se = NaN;
        ci = [t_last, Inf];
    end
end

function [se, ci, cross_rate] = bootstrap_pairs_timewin(t, r, w, q, sigma_base, kSigma, B, alpha, ...
                                                        tlo, thi, base_sec, plat_sec, baseFrac, platFrac)
    n = numel(t);
    t_last = t(end);
    tq = nan(B,1); crossed = false(B,1);

    for b = 1:B
        idx = randi(n, n, 1);
        t_b = t(idx); r_b = r(idx); w_b = w(idx);
        tlo_b = []; thi_b = [];
        if ~isempty(tlo) && ~isempty(thi)
            tlo_b = tlo(idx); thi_b = thi(idx);
        end

        [t_b, ord] = sort(t_b);
        r_b = r_b(ord); w_b = w_b(ord);
        if ~isempty(tlo_b), tlo_b = tlo_b(ord); end
        if ~isempty(thi_b), thi_b = thi_b(ord); end

        yhat_b = weighted_pav(r_b, w_b);

        widths_b = bin_widths(t_b, tlo_b, thi_b);
        [idxB, idxP] = select_time_windows(t_b, widths_b, base_sec, plat_sec, baseFrac, platFrac);

        b_b = wmedian(yhat_b(idxB), w_b(idxB));
        p_b = wmedian(yhat_b(idxP), w_b(idxP));
        A_b = p_b - b_b;
        theta_b = b_b + q*A_b;

        if A_b < kSigma * sigma_base
            tq(b) = t_last; crossed(b) = false;
        else
            [tb, cens] = crossing_time_from_fit(t_b, yhat_b, theta_b);
            tq(b) = tb; crossed(b) = ~cens;
        end
    end

    cross_rate = mean(crossed);
    if any(crossed)
        x = tq(crossed);
        se = std(x);
        lo = quantile(x, alpha/2);
        hi = quantile(x, 1 - alpha/2);
        ci = [lo, hi];
    else
        se = NaN;
        ci = [t_last, Inf];
    end
end

% ----------------- Window selection by TIME -----------------

function [idxB, idxP] = select_time_windows(t, widths, base_sec, plat_sec, baseFrac, platFrac)
% Choose baseline and plateau index sets by a target TIME length (seconds).
% If base_sec/plat_sec are empty, fall back to earliest/latest fractions.
    n = numel(t);
    if isempty(widths), widths = ones(n,1); end

    tmin = t(1); tmax = t(end);
    totalT = (tmax - tmin);

    if isempty(base_sec)
        kB = max(1, ceil(baseFrac * n));
        idxB = 1:kB;
    else
        cum = cumsum(widths);
        idxB = find(cum <= max(base_sec, 0));
        if isempty(idxB), idxB = 1; end
    end

    if isempty(plat_sec)
        kP = max(1, ceil(platFrac * n));
        idxP = (n-kP+1):n;
    else
        cum_rev = cumsum(flipud(widths));
        idx_rev = find(cum_rev <= max(plat_sec, 0));
        if isempty(idx_rev), idxP = n; else, idxP = n - idx_rev + 1; end
    end

    idxB = unique(idxB(:))';
    idxP = unique(idxP(:))';
end

function widths = bin_widths(t, tlo, thi)
% Estimate per-bin durations. Prefer edges; otherwise infer from centers.
    n = numel(t);
    widths = zeros(n,1);
    if ~isempty(tlo) && ~isempty(thi)
        widths = max(thi - tlo, 0);
        % fill any non-finite with center-based widths
        bad = ~isfinite(widths) | widths<=0;
        if any(bad)
            widths(bad) = center_widths(t(bad), t);
        end
    else
        widths = center_widths(t, t);
    end
    % ensure positive fallback
    widths(~isfinite(widths) | widths<=0) = eps;
end

function w = center_widths(tc, tall)
% Compute approximate bin widths from centers using neighbor midpoints.
    if nargin < 2, tall = tc; end
    t = tall(:);
    n = numel(t);
    % midpoints between neighbors
    mid = (t(1:end-1) + t(2:end))/2;
    left = [t(1) - (mid(1)-t(1)); mid];
    right = [mid; t(end) + (t(end)-mid(end))];
    % build a piecewise width look-up at centers in 't'
    base_w = right - left;

    % Map requested centers 'tc' onto the nearest position in 't'
    [~, loc] = ismember(tc(:), t);
    if any(loc==0)
        % if tc is not exactly in t, interpolate widths
        w = interp1(t, base_w, tc, 'nearest', 'extrap');
    else
        w = base_w(loc);
    end
end

% ----------------- Core helpers -----------------

function name = qfieldname(q)
    name = sprintf('q%g', round(100*q));
    name = strrep(name, '.', '_');
end

function opts = filldefaults(opts, defaults)
    f = fieldnames(defaults);
    for k = 1:numel(f)
        if ~isfield(opts, f{k}) || isempty(opts.(f{k}))
            opts.(f{k}) = defaults.(f{k});
        end
    end
end

function [t, r, w, tlo, thi] = extract_trw_edges(Si, useCountsW)
    have_t = isfield(Si,'t') && ~isempty(Si.t);
    have_edges = isfield(Si,'t_lo') && isfield(Si,'t_hi') ...
                 && ~isempty(Si.t_lo) && ~isempty(Si.t_hi);

    if have_t
        t = Si.t(:);
    elseif have_edges
        t = 0.5*(Si.t_lo(:) + Si.t_hi(:));
    else
        error('Sample is missing time information: provide t or (t_lo,t_hi).');
    end

    if ~isfield(Si,'r') || isempty(Si.r)
        error('Sample is missing response vector r.');
    end
    r = Si.r(:);

    if useCountsW && isfield(Si,'n') && ~isempty(Si.n)
        w = Si.n(:);
    else
        w = ones(size(r));
    end

    tlo = []; thi = [];
    if have_edges
        tlo = Si.t_lo(:); thi = Si.t_hi(:);
    end

    m = min([numel(t), numel(r), numel(w)]);
    t = t(1:m); r = r(1:m); w = w(1:m);
    if ~isempty(tlo), tlo = tlo(1:m); end
    if ~isempty(thi), thi = thi(1:m); end
end

function yhat = weighted_pav(y, w)
% Weighted Pool-Adjacent-Violators for nondecreasing fit.
    y = y(:); w = w(:);
    n = numel(y);
    mu   = y;
    ww   = w;
    idxs = num2cell(1:n);
    while true
        d = diff(mu);
        v = find(d < 0, 1, 'first');
        if isempty(v), break; end
        newW  = ww(v) + ww(v+1);
        newMu = (ww(v)*mu(v) + ww(v+1)*mu(v+1)) / newW;
        ww(v)   = newW;
        mu(v)   = newMu;
        ww(v+1) = [];
        mu(v+1) = [];
        idxs{v} = [idxs{v} idxs{v+1}];
        idxs(v+1) = [];
    end
    yhat = zeros(n,1);
    for b = 1:numel(mu)
        yhat(idxs{b}) = mu(b);
    end
end

function [tq, cens, idxCross] = crossing_time_from_fit(t, yhat, theta)
    j = find(yhat >= theta, 1, 'first');
    if isempty(j)
        tq = t(end); cens = true; idxCross = NaN;
    elseif j == 1
        tq = t(1);   cens = false; idxCross = 1;
    else
        t0 = t(j-1); t1 = t(j);
        y0 = yhat(j-1); y1 = yhat(j);
        if y1 <= y0
            tq = t1;
        else
            tq = t0 + (theta - y0) * (t1 - t0) / (y1 - y0);
        end
        cens = false; idxCross = j;
    end
end

function m = wmedian(x, w)
    x = x(:); w = w(:);
    [xs, idx] = sort(x);
    ws = w(idx); ws(~isfinite(ws)) = 0;
    tot = sum(ws);
    if tot <= 0, m = median(xs,'omitnan'); return; end
    c = cumsum(ws);
    k = find(c >= 0.5*tot, 1, 'first');
    m = xs(k);
end