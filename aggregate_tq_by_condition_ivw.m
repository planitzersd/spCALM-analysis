function T = aggregate_tq_by_condition_ivw(res, tagFields, alpha)
% aggregate_tq_by_condition_ivw
% Group replicate t_q estimates by condition tags and pool with
% inverse-variance weighting on the log-time scale (fixed-effect with
% over-dispersion inflation). Censored replicates are excluded.
%
% INPUTS
%   res       : struct array from tq_isotonic_struct, with fields at least:
%               - tq        (per-replicate t_q estimate; finite, >0)
%               - se        (bootstrap SE of t_q; may be NaN/empty)
%               - censored  (logical; true if replicate didn't reach threshold)
%               - reptag    (replicate id; not used for pooling but kept in notes)
%               - ... plus your condition tag fields, e.g. dose, temp, cellline, etc.
%   tagFields : cellstr of field names that define a condition (e.g., {'dose','cellline'})
%   alpha     : CI level (default 0.05 → 95% CI)
%
% OUTPUT
%   T : table, one row per unique condition, with columns:
%       [tag fields], k, n_censored_excluded, t_hat, se, ci_lo, ci_hi, ...
%       sum_w_log, Q, df, psi, method, note
%
% METHOD (per condition)
%   1) Keep replicates with censored==false, finite tq>0, finite se>0.
%   2) Work on z = log(tq), with se_z ≈ se / tq (delta method).
%   3) IVW pooled log-time: zhat = sum(w*z)/sum(w), w=1/se_z^2.
%   4) Over-dispersion inflation: Q = sum(w*(z-zhat)^2), df=k-1,
%      psi = max(1, Q/df); SEz = sqrt(psi / sum(w)).
%   5) CI on log scale with t_{df} (df>=1) else normal; back-transform.
%
% NOTES
%   - If k==0 → outputs NaNs with note 'all censored/invalid'.
%   - If k==1 → psi=1, uses normal critical value (z) for CI.
%   - If a tag field is missing for some rows, it is treated as "".
%
% Author: (you)

    if nargin < 3 || isempty(alpha), alpha = 0.05; end
    if ischar(tagFields), tagFields = {tagFields}; end

    N = numel(res);
    if N == 0
        T = table(); return
    end

    % Build condition keys and also capture displayed tag values
    keys = strings(N,1);
    tagVals = cell(numel(tagFields), 1);
    for j = 1:numel(tagFields), tagVals{j} = strings(N,1); end

    for i = 1:N
        parts = strings(1, numel(tagFields));
        for j = 1:numel(tagFields)
            fj = tagFields{j};
            if isfield(res, fj) && ~isempty(res(i).(fj))
                parts(j) = stringify_tagval(res(i).(fj));
                tagVals{j}(i) = parts(j);
            else
                parts(j) = "";
                tagVals{j}(i) = "";
            end
        end
        keys(i) = strjoin(parts, " | ");
    end

    [uKeys, ~, groupIdx] = unique(keys, 'stable');

    % Prepare output containers
    out = struct();
    outRows = numel(uKeys);
    for r = 1:outRows
        out(r).k = NaN; %#ok<*AGROW>
    end

    % Loop over conditions
    for g = 1:numel(uKeys)
        idx = find(groupIdx == g);
        % Extract replicate-level values and exclude censored/invalid
        tq    = getfield_default(res, idx, 'tq');
        se    = getfield_default(res, idx, 'se');
        cens  = getfield_default(res, idx, 'censored');

        valid = ~logical(cens) & isfinite(tq) & (tq > 0) & isfinite(se) & (se > 0);
        t_use = tq(valid);
        s_use = se(valid);

        out(g).k = numel(t_use);
        out(g).n_censored_excluded = sum(~valid);

        % Copy tag values (first non-empty per field)
        for j = 1:numel(tagFields)
            vals_g = tagVals{j}(idx);
            nz = find(vals_g ~= "", 1, 'first');
            if isempty(nz), v = ""; else, v = vals_g(nz); end
            out(g).(tagFields{j}) = v;
        end

        if out(g).k == 0
            % No usable replicates
            out(g).t_hat = NaN;
            out(g).se = NaN;
            out(g).ci_lo = NaN;
            out(g).ci_hi = NaN;
            out(g).sum_w_log = NaN;
            out(g).Q = NaN; out(g).df = 0; out(g).psi = NaN;
            out(g).method = "IVW log-time (no usable reps)";
            out(g).note = "all censored/invalid";
            continue
        end

        % Log-time work quantities
        z     = log(t_use);
        se_z  = s_use ./ t_use;            % delta method
        w     = 1 ./ (se_z.^2);
        sumw  = sum(w);
        zhat  = sum(w .* z) / sumw;

        % Over-dispersion
        if out(g).k >= 2
            Q  = sum(w .* (z - zhat).^2);
            df = out(g).k - 1;
            psi = max(1, Q / df);
            crit = tinv(1 - alpha/2, df);
        else
            Q = 0; df = 0; psi = 1;
            crit = norminv(1 - alpha/2);
        end
        SEz = sqrt(psi / sumw);

        % Back-transform
        t_hat = exp(zhat);
        ci_log = [zhat - crit*SEz, zhat + crit*SEz];
        ci = exp(ci_log);
        se_hat = t_hat * SEz;   % delta-method SE on original scale

        % Store
        out(g).t_hat = t_hat;
        out(g).se = se_hat;
        out(g).ci_lo = ci(1);
        out(g).ci_hi = ci(2);
        out(g).sum_w_log = sumw;
        out(g).Q = Q; out(g).df = df; out(g).psi = psi;
        out(g).method = "IVW log-time";
        out(g).note = "";
    end

    % Convert to table in a tidy order: tags first, then stats
    T = struct2table(out);
    % Reorder columns: tags, then counts/stats
    tagCols = tagFields(:)';
    statCols = {'k','n_censored_excluded','t_hat','se','ci_lo','ci_hi', ...
                'sum_w_log','Q','df','psi','method','note'};
    T = T(:, [tagCols, statCols]);
end

% ---------- helpers ----------
function s = stringify_tagval(v)
    % Convert various tag value types to a readable string cell scalar
    if isstring(v)
        s = v(1);
    elseif ischar(v)
        s = string(v);
    elseif islogical(v) && isscalar(v)
        s = string(mat2str(v));
    elseif isnumeric(v) && isscalar(v) && isfinite(v)
        s = string(num2str(v));
    else
        % fallback: JSON-like
        try
            s = string(jsonencode(v));
        catch
            s = "<complex>";
        end
    end
end

function vals = getfield_default(res, idx, fname)
    % Get field values as a column vector; default NaN/false if missing
    if ~isfield(res, fname)
        % choose default type by name
        if strcmp(fname,'censored'), vals = false(numel(idx),1);
        else, vals = NaN(numel(idx),1);
        end
        return
    end
    vals = arrayfun(@(k) res(k).(fname), idx, 'UniformOutput', false);
    try
        vals = vertcat(vals{:});
    catch
        % non-uniform -> coerce
        vals = cellfun(@(x) double(x), vals, 'UniformOutput', true);
    end
end
