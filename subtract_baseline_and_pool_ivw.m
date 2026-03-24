function [T_cond, T_rep] = subtract_baseline_and_pool_ivw(res_tq, res_t0, tagFields, alpha)
% Join by 'name', compute d = tq - t0 with independent SE propagation
% (se_d = sqrt(se_q^2 + t0_se^2)), drop censored tq, then pool within
% conditions using IVW (log-time if all d>0, else original).
%
% res_t0 fields:  name, t0, t0_se

    if nargin < 4 || isempty(alpha), alpha = 0.05; end
    if ischar(tagFields), tagFields = {tagFields}; end

    % ---- lookup: t0 by name (uses t0 and t0_se) ----
    t0map = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(res_t0)
        if ~isfield(res_t0(i),'name') || isempty(res_t0(i).name), continue; end
        nm   = char(res_t0(i).name);
        t0   = getnum(res_t0(i), 't0');
        se0  = getnum(res_t0(i), 't0_se');
        t0map(nm) = struct('t0', t0, 'se0', se0);
    end

    % ---- join, subtract, propagate SE (independent) ----
    rep = []; k = 0;
    for i = 1:numel(res_tq)
        if ~isfield(res_tq(i),'name') || isempty(res_tq(i).name), continue; end
        nm = char(res_tq(i).name);
        if ~isKey(t0map, nm), continue; end

        % tq, se_q, censored (allow a few common aliases for tq/se/censored)
        tq  = picknum(res_tq(i), {'tq','t_hat','t50','t_q'});
        seq = picknum(res_tq(i), {'se','se_tq'});
        cens = picklog(res_tq(i), {'censored','is_censored'}, true);

        % start row with defaults so table always has all variables
        k = k+1;
        rep(k).name   = string(nm); %#ok<AGROW>
        rep(k).reptag = picknum(res_tq(i), {'reptag'});

        % condition tag values
        for j = 1:numel(tagFields)
            fj = tagFields{j};
            rep(k).(fj) = val2str(getfield_def(res_tq(i), fj, ""));
        end

        % raw fields
        rep(k).tq   = tq;
        rep(k).se_q = seq;
        rep(k).t0   = t0map(nm).t0;
        rep(k).se_0 = t0map(nm).se0;

        % adjusted (predeclare)
        rep(k).d    = NaN;
        rep(k).se_d = NaN;

        % flags
        rep(k).used   = false;
        rep(k).reason = "init";

        % exclusions
        if cens
            rep(k).reason = "censored_tq"; continue
        end
        if ~isfinite(tq) || ~isfinite(seq) || seq<=0
            rep(k).reason = "invalid_tq_or_se"; continue
        end
        if ~isfinite(rep(k).t0) || ~isfinite(rep(k).se_0) || rep(k).se_0<=0
            rep(k).reason = "invalid_t0_or_t0_se"; continue
        end

        % adjusted value and SE (independence)
        rep(k).d    = tq - rep(k).t0;
        rep(k).se_d = sqrt(seq.^2 + rep(k).se_0.^2);
        rep(k).used   = true;
        rep(k).reason = "ok";
    end

    if isempty(rep)
        T_rep = table(); T_cond = table(); return
    end
    T_rep = struct2table(rep);

    % ---- group by condition tags ----
    keys = build_keys(T_rep, tagFields);
    [uKeys, ~, gIdx] = unique(keys, 'stable');

    % collect rows in a cell array to avoid dissimilar-struct errors
    rowsC = {};
    for g = 1:numel(uKeys)
        sel = (gIdx == g);

        row = struct();
        % tags
        for j = 1:numel(tagFields)
            fj = tagFields{j};
            vals = T_rep.(fj)(sel);
            nz = find(vals ~= "", 1, 'first');
            row.(fj) = iif(~isempty(nz), vals(nz), "");
        end

        % usable reps
        used  = T_rep.used(sel);
        d_use = T_rep.d(sel);
        se_use= T_rep.se_d(sel);

        row.k          = sum(used);
        row.n_excluded = sum(sel) - row.k;

        if row.k == 0
            % empty condition
            [row.method,row.t_hat,row.se,row.ci_lo,row.ci_hi,row.sum_w,row.Q,row.df,row.psi,row.note] = ...
                deal("none",NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,"no usable replicates");
            rowsC{end+1} = row; %#ok<AGROW>
            continue
        end

        d_use = d_use(used);
        se_use= se_use(used);

        if all(d_use > 0)
            % ---------- IVW on log-time ----------
            z    = log(d_use);
            se_z = se_use ./ d_use;
            w    = 1 ./ (se_z.^2);
            sumw = sum(w);
            zhat = sum(w.*z) / sumw;

            if row.k >= 2
                Q  = sum(w.*(z - zhat).^2);
                df = row.k - 1;
                psi = max(1, Q/df);
                crit = tinv(1 - alpha/2, df);
            else
                Q = 0; df = 0; psi = 1;
                crit = norminv(1 - alpha/2);
            end
            SEz   = sqrt(psi / sumw);
            ciLog = [zhat - crit*SEz, zhat + crit*SEz];

            row.method = "IVW log-time";
            row.t_hat  = exp(zhat);
            row.se     = row.t_hat * SEz;
            row.ci_lo  = exp(ciLog(1));
            row.ci_hi  = exp(ciLog(2));
            row.sum_w  = sumw; row.Q = Q; row.df = df; row.psi = psi;
            row.note   = "";
        else
            % ---------- IVW on original scale ----------
            ok = isfinite(d_use) & isfinite(se_use) & se_use>0;
            d_use = d_use(ok); se_use = se_use(ok);
            kpos = numel(d_use);
            if kpos == 0
                [row.method,row.t_hat,row.se,row.ci_lo,row.ci_hi,row.sum_w,row.Q,row.df,row.psi,row.note] = ...
                    deal("none",NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,"no usable (nonnegative) d");
                rowsC{end+1} = row; %#ok<AGROW>
                continue
            end
            w    = 1 ./ (se_use.^2);
            sumw = sum(w);
            that = sum(w.*d_use) / sumw;

            if kpos >= 2
                Q  = sum(w.*(d_use - that).^2);
                df = kpos - 1;
                psi = max(1, Q/df);
                crit = tinv(1 - alpha/2, df);
            else
                Q = 0; df = 0; psi = 1;
                crit = norminv(1 - alpha/2);
            end
            SEt = sqrt(psi / sumw);

            row.method = "IVW original";
            row.t_hat  = that;
            row.se     = SEt;
            row.ci_lo  = that - crit*SEt;
            row.ci_hi  = that + crit*SEt;
            row.sum_w  = sumw; row.Q = Q; row.df = df; row.psi = psi;
            row.note   = "fallback (some d<=0)";
        end

        rowsC{end+1} = row; %#ok<AGROW>
    end

    % build condition table
    if isempty(rowsC)
        T_cond = table();
    else
        rows = [rowsC{:}];                 % now a uniform struct array
        T_cond = struct2table(rows);
        T_cond = T_cond(:, [tagFields(:)' {'k','n_excluded','method','t_hat','se','ci_lo','ci_hi','sum_w','Q','df','psi','note'}]);
    end
end

% ================= helpers =================
function x = getfield_def(S, f, def)
    if isfield(S,f) && ~isempty(S.(f)), x = S.(f); else, x = def; end
end
function x = val2str(v)
    if isstring(v), x = v(1);
    elseif ischar(v), x = string(v);
    elseif isnumeric(v) && isscalar(v) && isfinite(v), x = string(num2str(v));
    elseif islogical(v) && isscalar(v), x = string(mat2str(v));
    else, x = "";
    end
end
function x = picknum(S, cand)
    x = NaN;
    for i = 1:numel(cand)
        f = cand{i};
        if isfield(S,f) && ~isempty(S.(f)) && isnumeric(S.(f))
            x = S.(f); return
        end
    end
end
function x = picklog(S, cand, def)
    x = def;
    for i = 1:numel(cand)
        f = cand{i};
        if isfield(S,f) && ~isempty(S.(f))
            x = logical(S.(f)); return
        end
    end
end
function x = getnum(S, f)
    if isfield(S,f) && ~isempty(S.(f)) && isnumeric(S.(f))
        x = S.(f);
    else
        x = NaN;
    end
end
function keys = build_keys(T, tagFields)
    N = height(T);
    parts = strings(N, numel(tagFields));
    for j = 1:numel(tagFields)
        v = T.(tagFields{j});
        if ~isstring(v), v = string(v); end
        parts(:,j) = v;
    end
    keys = strings(N,1);
    for i = 1:N, keys(i) = strjoin(parts(i,:), " | "); end
end
function y = iif(cond, a, b)
    if cond, y = a; else, y = b; end
end
