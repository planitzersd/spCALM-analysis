function [T_samples, T_conditions] = summarize_delta_t50_t0(OUT, tag_fields, target_folder, opts)
% SUMMARIZE_DELTA_T50_T0
% Compute per-sample and per-condition statistics for Δ = t50(4PL) - t0(event-rate),
% with uncertainty propagation and inverse-variance weighted condition means.

if nargin < 4, opts = struct; end
if nargin < 3, target_folder = ''; end
if ischar(tag_fields) || isstring(tag_fields), tag_fields = cellstr(tag_fields); end

% ---- defaults
d.alpha       = 0.05;
d.t0_se_mode  = 'auto';       % 'auto'|'binwidth'|'fixed'|'zero'
d.t0_se_value = 0;
if isfield(OUT,'opts') && isfield(OUT.opts,'t0_binwidth')
    d.binwidth = OUT.opts.t0_binwidth;
else
    d.binwidth = 1.0;
end
fn = fieldnames(d);
for i=1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end
end
z = norminv(1 - opts.alpha/2);

S = OUT.samples;
Ns = numel(S);

% ---- per-sample extraction
t50   = nan(Ns,1);
se50  = nan(Ns,1);
t0    = nan(Ns,1);
se0   = nan(Ns,1);
has4  = false(Ns,1);

for s = 1:Ns
    % t0 and its SE
    if isfield(S(s),'t0') && ~isempty(S(s).t0) && isfinite(S(s).t0)
        t0(s) = S(s).t0;
        % preferred: sample-provided SE
        if strcmpi(opts.t0_se_mode, 'auto')
            if isfield(S(s),'t0_se') && ~isempty(S(s).t0_se) && isfinite(S(s).t0_se)
                se0(s) = max(0, S(s).t0_se);
            else
                se0(s) = opts.binwidth / sqrt(12);
            end
        else
            switch lower(opts.t0_se_mode)
                case 'zero'
                    se0(s) = 0;
                case 'fixed'
                    se0(s) = max(0, opts.t0_se_value);
                case 'binwidth'
                    se0(s) = opts.binwidth / sqrt(12);
                otherwise
                    se0(s) = opts.binwidth / sqrt(12);
            end
        end
    else
        t0(s)  = NaN;
        se0(s) = NaN;
    end

    % 4PL t50 and SE
    if isfield(S(s),'logistic4') && ~isempty(S(s).logistic4) ...
            && isfield(S(s).logistic4,'params') && isfield(S(s).logistic4.params,'t50')
        t50(s) = S(s).logistic4.params.t50;
        if isfield(S(s).logistic4,'se') && isfield(S(s).logistic4.se,'t50') && ~isempty(S(s).logistic4.se.t50)
            se50(s) = S(s).logistic4.se.t50;
        elseif isfield(S(s).logistic4,'ci') && isfield(S(s).logistic4.ci,'t50') && ~isempty(S(s).logistic4.ci.t50)
            ci = S(s).logistic4.ci.t50(:);
            if numel(ci) >= 2 && all(isfinite(ci))
                se50(s) = (ci(2) - ci(1)) / (2*z);
            end
        end
        has4(s) = true;
    end
end

Delta  = t50 - t0;
VarDel = se50.^2 + se0.^2;
SEdel  = sqrt(VarDel);
CIlo   = Delta - z*SEdel;
CIhi   = Delta + z*SEdel;

% ---- tags table
Tags = table();
for i = 1:numel(tag_fields)
    f = tag_fields{i};
    vals = strings(Ns,1);
    for s = 1:Ns
        if isfield(S(s), f) && ~isempty(S(s).(f)), vals(s) = string(S(s).(f)); else, vals(s) = "MISSING"; end
    end
    Tags.(f) = vals;
end
Name = strings(Ns,1);
for s=1:Ns
    if isfield(S(s),'name') && ~isempty(S(s).name), Name(s) = string(S(s).name); else, Name(s) = "sample_"+s; end
end

T_samples = [ table(Name, t50, se50, t0, se0, Delta, SEdel, CIlo, CIhi), Tags ];

% ---- per-condition inverse-variance weighted mean (exclude reptag)
key_fields = tag_fields(~strcmpi(tag_fields,'reptag'));
if isempty(key_fields)
    cond_key = repmat("ALL", Ns, 1);
else
    cond_key = join(Tags{:, key_fields}, "__", 2);
end
cond_key(ismissing(cond_key)) = "MISSING";
[uniqK, ~, Kidx] = unique(cond_key);
Nc = numel(uniqK);

CondRows = cell(Nc,1);
for k = 1:Nc
    rows = find(Kidx == k & has4 & isfinite(Delta) & isfinite(SEdel) & SEdel > 0);
    baseRow = table();
    for i = 1:numel(key_fields)
        f = key_fields{i};
        baseRow.(f) = Tags.(f)( max(rows(1),1) );
    end
    if isempty(rows)
        Row = [ baseRow, table(0, NaN, NaN, NaN, NaN, 'VariableNames', {'Nrep','meanDelta','SE','CIlo','CIhi'}) ];
    else
        w = 1 ./ (SEdel(rows).^2);
        w(~isfinite(w)) = 0;
        if sum(w) == 0
            m = mean(Delta(rows),'omitnan');
            se = std(Delta(rows),'omitnan') / sqrt(numel(rows));
        else
            m = sum(w .* Delta(rows)) / sum(w);
            se = sqrt(1 / sum(w));
        end
    Row = [ baseRow, table(numel(rows), m, se, m - z*se, m + z*se, 'VariableNames', {'Nrep','meanDelta','SE','CIlo','CIhi'}) ];
    end
    CondRows{k} = Row;
end
if Nc > 0
    T_conditions = vertcat(CondRows{:});
    if ~isempty(key_fields)
        T_conditions = sortrows(T_conditions, key_fields);
    end
else
    T_conditions = table();
end

if ~isempty(target_folder)
    if ~exist(target_folder,'dir'), mkdir(target_folder); end
    writetable(T_samples,    fullfile(target_folder, 'delta_t50_minus_t0__samples.csv'));
    writetable(T_conditions, fullfile(target_folder, 'delta_t50_minus_t0__conditions.csv'));
end
end
