function [Tlog_samp, Tlog_cond, Tgam_samp, Tgam_cond, Tbest_samp, Tbest_cond] = ...
    summarize_fit_parameters(OUT, tag_fields, target_folder, opts)
% SUMMARIZE_FIT_PARAMETERS
% Build per-sample and per-condition summaries (with uncertainty) for:
%   - Logistic 4PL: B, A, top=B+A, span=A, t50, k
%   - Gamma onset : B, A, top=B+A, span=A, t0, N, k
%   - Best span   : choose span (and SE) from lower BIC (logistic vs gamma)
%
% Inputs
%   OUT          : struct from fit_binned_kinetics (OUT.samples(:).logistic4 / .gamma / .best)
%   tag_fields   : cellstr of fields defining conditions (e.g., {'pHtag','NAItag','recsattag','reptag'})
%   target_folder: (optional) folder to save CSVs (''/[] => no writing)
%   opts         : (optional) struct:
%                  .alpha (default 0.05) for 95% CI
%
% Outputs
%   Tlog_samp, Tlog_cond : logistic per-sample and per-condition tables
%   Tgam_samp, Tgam_cond : gamma    per-sample and per-condition tables
%   Tbest_samp,Tbest_cond: best-span per-sample and per-condition tables
%
% Notes
% - If SE is missing but 95% CI exists, SE is approximated as (CI_hi - CI_lo)/(2*z).
% - For 'top' we assume Cov(B,A) ~ 0, so SE_top^2 = SE_B^2 + SE_A^2.
% - Per-condition stats: inverse-variance weighted mean; fallback to unweighted if SEs missing.

    if nargin < 4, opts = struct; end
    if nargin < 3, target_folder = ''; end
    if ischar(tag_fields) || isstring(tag_fields), tag_fields = cellstr(tag_fields); end

    % -------- settings --------
    d.alpha = 0.05;
    fn = fieldnames(d);
    for i=1:numel(fn), if ~isfield(opts,fn{i}), opts.(fn{i}) = d.(fn{i}); end, end
    z = norminv(1 - opts.alpha/2);

    S = OUT.samples;
    Ns = numel(S);

    % -------- helpers --------
    function se = se_from_ci(ci)
        se = NaN;
        if ~isempty(ci)
            v = ci(:);
            if numel(v) >= 2 && all(isfinite(v(1:2)))
                se = (v(2) - v(1)) / (2*z);
            end
        end
    end

    function [vals, SEs, ok] = pull_param(model, pname)
        vals = nan(Ns,1); SEs = nan(Ns,1); ok = false(Ns,1);
        for s=1:Ns
            if ~isfield(S(s), model) || isempty(S(s).(model)), continue; end
            M = S(s).(model);
            if isfield(M,'params') && isfield(M.params, pname)
                vals(s) = M.params.(pname);
                if isfield(M,'se') && isfield(M.se, pname) && ~isempty(M.se.(pname))
                    SEs(s) = M.se.(pname);
                elseif isfield(M,'ci') && isfield(M.ci, pname) && ~isempty(M.ci.(pname))
                    SEs(s) = se_from_ci(M.ci.(pname));
                end
                ok(s) = isfinite(vals(s));
            end
        end
    end

    function bic = pull_bic(model)
        bic = inf(Ns,1);
        for s=1:Ns
            if isfield(S(s), model) && ~isempty(S(s).(model)) && isfield(S(s).(model),'bic')
                b = S(s).(model).bic;
                if isfinite(b), bic(s) = b; end
            end
        end
    end

    function Ttags = build_tags()
        Ttags = table();
        for i = 1:numel(tag_fields)
            f = tag_fields{i};
            vals = strings(Ns,1);
            for s = 1:Ns
                if isfield(S(s), f) && ~isempty(S(s).(f)), vals(s) = string(S(s).(f)); else, vals(s) = "MISSING"; end
            end
            Ttags.(f) = vals;
        end
    end

    function Name = get_names()
        Name = strings(Ns,1);
        for s=1:Ns
            if isfield(S(s),'name') && ~isempty(S(s).name), Name(s) = string(S(s).name);
            else, Name(s) = "sample_"+s; end
        end
    end

    function [Tsamp, Tcond] = make_tables(param_names, model_name, Ttags)
        % param_names: cellstr of parameter base names and derived 'top','span'

        % extract parameters & SEs
        data = struct();
        for i=1:numel(param_names)
            p = param_names{i};
            if strcmpi(p,'top') || strcmpi(p,'span')
                % will compute later
                continue
            end
            [vals, SEs, ok] = pull_param(model_name, p);
            data.(p)    = vals;
            data.([p '_se']) = SEs;
            data.([p '_ok']) = ok;
        end

        % derived: top = B + A ; span = A
        if isfield(data,'B') && isfield(data,'A')
            top = data.B + data.A;
            top_se = sqrt( nansqr(data.B_se) + nansqr(data.A_se) ); % independence assumption
            data.top = top; data.top_se = top_se;
        end
        if isfield(data,'A')
            data.span = data.A; data.span_se = data.A_se;
        end

        % per-sample table
        Name = get_names();
        Tsamp = table(Name);
        base_order = {}; % keep output order nice
        for i=1:numel(param_names)
            p = param_names{i};
            v = nan(Ns,1); se = nan(Ns,1);
            if isfield(data, p), v = data.(p); end
            if isfield(data, [p '_se']), se = data.([p '_se']); end
            Tsamp.(p)     = v;
            Tsamp.([p '_SE'])  = se;
            Tsamp.([p '_CIlo'])= v - z.*se;
            Tsamp.([p '_CIhi'])= v + z.*se;
            base_order = [base_order, {p, [p '_SE'], [p '_CIlo'], [p '_CIhi']}]; %#ok<AGROW>
        end
        Tsamp = [Tsamp Ttags];

        % per-condition aggregation (exclude reptag)
        key_fields = tag_fields(~strcmpi(tag_fields,'reptag'));
        if isempty(key_fields)
            cond_key = repmat("ALL", Ns, 1);
        else
            cond_key = join(Ttags{:, key_fields}, "__", 2);
        end
        cond_key(ismissing(cond_key)) = "MISSING";
        [uniqK, ~, Kidx] = unique(cond_key);
        Nc = numel(uniqK);

        % build rows param-by-param with inverse-variance weighted mean
        CondRows = cell(Nc,1);
        for k = 1:Nc
            rows = find(Kidx == k);
            baseRow = table();
            for i = 1:numel(key_fields)
                f = key_fields{i};
                baseRow.(f) = Ttags.(f)( max(rows(1),1) );
            end

            % For each param, compute mean/SE/CI across rows
            Row = baseRow;
            for i=1:numel(param_names)
                p = param_names{i};
                v = Tsamp.(p)(rows);
                se = Tsamp.([p '_SE'])(rows);

                use = isfinite(v) & isfinite(se) & se > 0;
                if any(use)
                    w = 1 ./ (se(use).^2);
                    m = sum(w .* v(use)) / sum(w);
                    se_m = sqrt( 1 / sum(w) );
                else
                    % fallback if no SEs: unweighted mean & SE by sample SD
                    v2 = v(isfinite(v));
                    if ~isempty(v2)
                        m = mean(v2);
                        se_m = std(v2) / sqrt(numel(v2));
                    else
                        m = NaN; se_m = NaN;
                    end
                end
                Row.([p '_mean'])  = m;
                Row.([p '_SE'])    = se_m;
                Row.([p '_CIlo'])  = m - z*se_m;
                Row.([p '_CIhi'])  = m + z*se_m;
            end
            % Also include Nrep used per param? For simplicity, a single Nrep (rows with finite p)
            Nrep = sum(isfinite(Tsamp.(param_names{1})(rows)));
            Row.Nrep = Nrep;

            CondRows{k} = Row;
        end
        if Nc > 0
            Tcond = vertcat(CondRows{:});
            if ~isempty(key_fields)
                Tcond = sortrows(Tcond, key_fields);
            end
        else
            Tcond = table();
        end
    end

    function y = nansqr(x), y = x.^2; y(~isfinite(y)) = NaN; end

    % -------- build tags once --------
    Ttags = build_tags();

    % -------- LOGISTIC --------
    log_params = {'B','A','top','span','t50','k'};
    [Tlog_samp, Tlog_cond] = make_tables(log_params, 'logistic4', Ttags);

    % -------- GAMMA --------
    gam_params = {'B','A','top','span','t0','N','k'};
    [Tgam_samp, Tgam_cond] = make_tables(gam_params, 'gamma', Ttags);

    % -------- BEST span per sample (by lower BIC) --------
    bicL = pull_bic('logistic4');
    bicG = pull_bic('gamma');

    chooseLog = bicL <= bicG;
    % Per-sample best span and SE
    spanL = Tlog_samp.span;      seL = Tlog_samp.span_SE;
    spanG = Tgam_samp.span;      seG = Tgam_samp.span_SE;

    bestSpan  = nan(Ns,1);
    bestSE    = nan(Ns,1);
    bestModel = strings(Ns,1);
    for s=1:Ns
        if isfinite(bicL(s)) || isfinite(bicG(s))
            if chooseLog(s)
                bestSpan(s) = spanL(s); bestSE(s) = seL(s); bestModel(s) = "logistic4";
            else
                bestSpan(s) = spanG(s); bestSE(s) = seG(s); bestModel(s) = "gamma";
            end
        else
            bestModel(s) = "MISSING";
        end
    end
    Name = get_names();
    Tbest_samp = table(Name, bestModel, bestSpan, bestSE, ...
                       bestSpan - z*bestSE, bestSpan + z*bestSE, ...
                       'VariableNames', {'Name','bestModel','span','span_SE','span_CIlo','span_CIhi'});
    Tbest_samp = [Tbest_samp Ttags];

    % Per-condition combine best span
    key_fields = tag_fields(~strcmpi(tag_fields,'reptag'));
    if isempty(key_fields)
        cond_key = repmat("ALL", Ns, 1);
    else
        cond_key = join(Ttags{:, key_fields}, "__", 2);
    end
    cond_key(ismissing(cond_key)) = "MISSING";
    [uniqK, ~, Kidx] = unique(cond_key);
    Nc = numel(uniqK);

    CondRows = cell(Nc,1);
    for k = 1:Nc
        rows = find(Kidx == k);
        baseRow = table();
        for i = 1:numel(key_fields)
            f = key_fields{i};
            baseRow.(f) = Ttags.(f)( max(rows(1),1) );
        end
        v  = Tbest_samp.span(rows);
        se = Tbest_samp.span_SE(rows);

        use = isfinite(v) & isfinite(se) & se > 0;
        if any(use)
            w = 1 ./ (se(use).^2);
            m = sum(w .* v(use)) / sum(w);
            se_m = sqrt( 1 / sum(w) );
        else
            v2 = v(isfinite(v));
            if ~isempty(v2)
                m = mean(v2);
                se_m = std(v2) / sqrt(numel(v2));
            else
                m = NaN; se_m = NaN;
            end
        end
        Row = baseRow;
        Row.Nrep      = sum(isfinite(v));
        Row.span_mean = m;
        Row.span_SE   = se_m;
        Row.span_CIlo = m - z*se_m;
        Row.span_CIhi = m + z*se_m;
        CondRows{k} = Row;
    end
    if Nc > 0
        Tbest_cond = vertcat(CondRows{:});
        if ~isempty(key_fields)
            Tbest_cond = sortrows(Tbest_cond, key_fields);
        end
    else
        Tbest_cond = table();
    end

    % -------- write CSVs if requested --------
    if ~isempty(target_folder)
        if ~exist(target_folder,'dir'), mkdir(target_folder); end
        writetable(Tlog_samp,  fullfile(target_folder, 'logistic_params__samples.csv'));
        writetable(Tlog_cond,  fullfile(target_folder, 'logistic_params__conditions.csv'));
        writetable(Tgam_samp,  fullfile(target_folder, 'gamma_params__samples.csv'));
        writetable(Tgam_cond,  fullfile(target_folder, 'gamma_params__conditions.csv'));
        writetable(Tbest_samp, fullfile(target_folder, 'best_span__samples.csv'));
        writetable(Tbest_cond, fullfile(target_folder, 'best_span__conditions.csv'));
    end
end
