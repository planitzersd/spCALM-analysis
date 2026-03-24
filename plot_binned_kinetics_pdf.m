function plot_binned_kinetics_pdf(OUT, pdf_file, opts)
% PLOT_BINNED_KINETICS_PDF2
% Letter-sized multi-page PDF with 3x3 samples per page.
% Each sample uses two tiles: left (data + BOTH model curves, solid lines),
% right (compact parameter panel with 3 sig figs, CIs, R^2/AIC/BIC).
%
% Usage:
%   plot_binned_kinetics_pdf2(OUT, 'report.pdf');
%   plot_binned_kinetics_pdf2(OUT, 'report.pdf', struct('shade_pre_t0',true));
%
% INPUTS:
%   OUT      : struct from your fit pipeline (expects samples(:).bins, .logistic4, .gamma, .best)
%   pdf_file : output PDF path (overwrites first, then appends)
%   opts     : optional struct:
%       .shade_pre_t0      (default true)      shade t<t0 if samples(s).t0 exists
%       .points_per_curve  (default 400)       curve resolution
%       .base_font         (default 9)         font size for axes
%       .panel_font        (default 8)         font size for panel text
%       .errorbar_color    (default [0.2 0.2 0.2])
%       .log_color         (default [0 0 0])        % logistic 4PL
%       .gam_color         (default [0 0.40 0.80])  % gamma onset
%       .shade_color       (default [0.93 0.95 1.00])
%       .panel_bg          (default [1 1 1])
%       .panel_edge        (default [0.85 0.85 0.85])
%
% OUTPUT:
%   Writes/append pages to pdf_file.

% -------- defaults --------
if nargin < 3, opts = struct; end
D = struct( ...
    'shade_pre_t0', true, ...
    'points_per_curve', 400, ...
    'base_font', 9, ...
    'panel_font', 8, ...
    'errorbar_color', [0.2 0.2 0.2], ...
    'log_color', [0 0 0], ...
    'gam_color', [0 0.40 0.80], ...
    'shade_color', [0.93 0.95 1.00], ...
    'panel_bg', [1 1 1], ...
    'panel_edge', [0.85 0.85 0.85] ...
);
fns = fieldnames(D);
for i=1:numel(fns)
    if ~isfield(opts, fns{i}), opts.(fns{i}) = D.(fns{i}); end
end

% -------- filter plottable samples --------
S = OUT.samples;
Ns = numel(S);
keep = false(Ns,1);
for s=1:Ns
    keep(s) = isfield(S(s),'bins') && ~isempty(S(s).bins) && ...
              (isfield(S(s),'logistic4') || isfield(S(s),'gamma'));
end
idx = find(keep);
if isempty(idx)
    warning('No samples with bins and at least one model.'); return;
end

% -------- page geometry: 3 rows x 6 columns (2 tiles per sample) --------
nrows = 3; ncols = 2; % => 3 samples per row, 9 per page
samples_per_page = nrows * (ncols/2); % 9
if exist(pdf_file,'file'), delete(pdf_file); end

% Letter size (inches)
pageW = 8.5; pageH = 11.0;

% -------- iterate pages --------
K = numel(idx);
pages = ceil(K / samples_per_page);
k0 = 1;

for pg = 1:pages
    k1 = min(K, k0 + samples_per_page - 1);
    idx_page = idx(k0:k1);

    % Figure: Letter-sized, Arial defaults
    f = figure('Color','w', 'Units','inches', ...
               'Position',[0.5 0.5 pageW pageH], ...
               'PaperUnits','inches', 'PaperPosition',[0 0 pageW pageH], ...
               'Renderer','painters');
    set(f, 'DefaultAxesFontName','Arial');
    set(f, 'DefaultTextFontName','Arial');

    tlo = tiledlayout(f, nrows, ncols, 'TileSpacing','compact', 'Padding','compact');
    tlo.TileIndexing = 'rowmajor';

    for kk = 1:numel(idx_page)
        s = idx_page(kk);

        % left tile: plot
        ax = nexttile(tlo, (kk-1)*2 + 1); % tiles: 1,3,5,... within row
        ax.FontSize = opts.base_font;
        hold(ax, 'on'); grid(ax, 'on');

        % Data table
        T = S(s).bins;
        if isempty(T) || height(T) < 1
            title(ax, get_sample_title(S, s), 'Interpreter','none');
            axis(ax,'off'); % also create panel but empty
        else
            t_min = min(T.t_lo); t_max = max(T.t_hi);
            tg = linspace(t_min, t_max, opts.points_per_curve).';

            % Optional pre-t0 shade
            if opts.shade_pre_t0 && isfield(S(s),'t0') && ~isempty(S(s).t0) && isfinite(S(s).t0)
                t0 = S(s).t0;
                xl = [t_min, min(t0,t_max), min(t0,t_max), t_min];
                yl = [0, 0, 1, 1];
                if xl(1) < xl(2)
                    patch('XData',xl,'YData',yl,'Parent',ax, ...
                          'FaceColor',opts.shade_color,'EdgeColor','none','FaceAlpha',0.35);
                    xline(ax, t0, ':', 'Color', opts.gam_color, 'LineWidth', 1.0);
                end
            end

            % Binned points with Wilson CI
            [lo, hi] = wilson_ci(T.p, T.n, 0.95);
            eb = errorbar(ax, T.t, T.p, T.p - lo, hi - T.p, 'o', ...
                          'MarkerSize', 3.5, 'LineWidth', 0.9, 'Color', opts.errorbar_color);
            eb.CapSize = 0;

            % Logistic 4PL (solid, black)
            if isfield(S(s),'logistic4') && ~isempty(S(s).logistic4)
                P = S(s).logistic4.params;
                fL = P.bottom + P.span ./ (1 + exp(-P.k .* (tg - P.t50)));
                plot(ax, tg, fL, '-', 'LineWidth', 1.6, 'Color', opts.log_color, 'DisplayName','4PL');
            end

            % Gamma-onset (solid, blue)
            if isfield(S(s),'gamma') && ~isempty(S(s).gamma)
                P = S(s).gamma.params;
                tau = max(tg - P.t0, 0);
                fG  = P.bottom + P.span .* gamcdf(tau, P.N, 1/max(P.k, eps));
                plot(ax, tg, fG, '-', 'LineWidth', 1.6, 'Color', opts.gam_color, 'DisplayName','Gamma');
            end

            ylim(ax,[0 1]); xlim(ax,[t_min t_max]);
            xlabel(ax,'Time'); ylabel(ax,'Fused fraction');
            title(ax, get_sample_title(S, s), 'Interpreter','none');

            % legend on first sample of the page only
            if kk==1
                legend(ax, 'Location','best');
            end
        end

        hold(ax, 'off');

        % right tile: parameter panel
        pax = nexttile(tlo, (kk-1)*2 + 2);
        pax.Visible = 'off'; % hide ticks; we'll only use rectangle + text
        % draw a light panel background using axes limits [0,1]
        set(pax, 'Visible','on', 'XLim',[0 1], 'YLim',[0 1], 'XColor','none', 'YColor','none');
        rectangle('Parent',pax, 'Position',[0 0 1 1], ...
                  'FaceColor',opts.panel_bg, 'EdgeColor',opts.panel_edge, 'LineWidth', 0.75);
        panel_str = build_panel_text(S(s));
        text(0.02, 0.98, panel_str, 'Parent',pax, 'Units','normalized', ...
             'HorizontalAlignment','left', 'VerticalAlignment','top', ...
             'FontName','Arial', 'FontSize', opts.panel_font, 'Color',[0 0 0]);
        pax.Visible = 'off'; % keep background/text only
    end

    % Export page
    try
        exportgraphics(f, pdf_file, 'ContentType','vector', 'Resolution',300, ...
                       'BackgroundColor','white', 'Append', (pg>1));
    catch
        print(f, '-dpdf', '-painters', pdf_file, '-bestfit', ['-r' num2str(300)]);
    end
    close(f);

    k0 = k1 + 1;
end
end

% ============================= helpers ==============================

function ttl = get_sample_title(S, s)
ttl = sprintf('Sample %d', s);
if isfield(S(s),'name') && ~isempty(S(s).name)
    ttl = string(S(s).name);
end
end

function panel_str = build_panel_text(Ss)
% 3 sig figs; include best model (by BIC), params, CIs, R^2/AIC/BIC
lines = strings(0,1);

% Which is best?
best = '';
if isfield(Ss,'best') && ~isempty(Ss.best)
    best = lower(string(Ss.best));
else
    best = choose_best(Ss);
end
if ~isempty(best)
    lines(end+1) = "Best (BIC): " + upper(strrep(best,'logistic4','4PL'));
end

% 4PL
if isfield(Ss,'logistic4') && ~isempty(Ss.logistic4)
    P = Ss.logistic4.params; C = getfield2(Ss.logistic4,'ci',[]);
    lines(end+1) = "4PL";
    lines(end+1) = sprintf('  B=%s  A=%s', g3(P.bottom), g3(P.span));
    lines(end+1) = sprintf('  t50=%s  k=%s', g3(P.t50), g3(P.k));
    if ~isempty(C)
        lines(end+1) = sprintf('  CI B=[%s,%s]  A=[%s,%s]', g3(C.B(1)), g3(C.B(2)), g3(C.A(1)), g3(C.A(2)));
        lines(end+1) = sprintf('  CI t50=[%s,%s]  k=[%s,%s]', g3(C.t50(1)), g3(C.t50(2)), g3(C.k(1)), g3(C.k(2)));
    end
    [R2, AIC, BIC] = quick_gof(Ss,'logistic4');
    lines(end+1) = sprintf('  R^2=%s  AIC=%s  BIC=%s', g3(R2), g3(AIC), g3(BIC));
end

% Gamma
if isfield(Ss,'gamma') && ~isempty(Ss.gamma)
    P = Ss.gamma.params; C = getfield2(Ss.gamma,'ci',[]);
    lines(end+1) = "Gamma";
    lines(end+1) = sprintf('  B=%s  A=%s', g3(P.bottom), g3(P.span));
    lines(end+1) = sprintf('  t0=%s  N=%s  k=%s', g3(P.t0), g3(P.N), g3(P.k));
    if ~isempty(C)
        lines(end+1) = sprintf('  CI B=[%s,%s]  A=[%s,%s]', g3(C.B(1)), g3(C.B(2)), g3(C.A(1)), g3(C.A(2)));
        lines(end+1) = sprintf('  CI t0=[%s,%s]  N=[%s,%s]  k=[%s,%s]', ...
            g3(C.t0(1)), g3(C.t0(2)), g3(C.N(1)), g3(C.N(2)), g3(C.k(1)), g3(C.k(2)));
    end
    [R2, AIC, BIC] = quick_gof(Ss,'gamma');
    lines(end+1) = sprintf('  R^2=%s  AIC=%s  BIC=%s', g3(R2), g3(AIC), g3(BIC));
end

panel_str = strjoin(lines, newline);
end

function best = choose_best(Ss)
hasL = isfield(Ss,'logistic4') && ~isempty(Ss.logistic4) && isfield(Ss.logistic4,'bic');
hasG = isfield(Ss,'gamma')     && ~isempty(Ss.gamma)     && isfield(Ss.gamma,'bic');
if hasL && hasG
    best = ternary(Ss.logistic4.bic <= Ss.gamma.bic, 'logistic4', 'gamma');
elseif hasL
    best = 'logistic4';
elseif hasG
    best = 'gamma';
else
    best = '';
end
best = lower(string(best));
end

function [R2, AIC, BIC] = quick_gof(Ss, modelname)
R2 = NaN; AIC = NaN; BIC = NaN;
if ~isfield(Ss,'bins') || isempty(Ss.bins), return; end
if ~isfield(Ss, modelname) || isempty(Ss.(modelname)) || ~isfield(Ss.(modelname),'fit'), return; end
T = Ss.bins;
if ~all(ismember({'p','n'}, T.Properties.VariableNames)) || ~isfield(Ss.(modelname).fit,'p_hat')
    return;
end
y  = T.p(:).*T.n(:);
ph = Ss.(modelname).fit.p_hat(:);
n  = T.n(:);
L  = min([numel(y), numel(n), numel(ph)]);
y  = y(1:L); n = n(1:L); ph = ph(1:L);
ll_model = sum( y .* log(max(ph,eps)) + (n - y) .* log(max(1 - ph, eps)) );
p_bar    = sum(y) / max(sum(n), eps);
ll_null  = sum( y .* log(max(p_bar,eps)) + (n - y) .* log(max(1 - p_bar,eps)) );
if isfinite(ll_model) && isfinite(ll_null) && ll_null ~= 0
    R2 = 1 - (ll_model / ll_null);
end
if isfield(Ss.(modelname),'aic'), AIC = Ss.(modelname).aic; end
if isfield(Ss.(modelname),'bic'), BIC = Ss.(modelname).bic; end
end

function [lo, hi] = wilson_ci(p, n, conf)
if nargin<3, conf=0.95; end
z = norminv(0.5+conf/2);
p = p(:); n = n(:);
n_eff = max(n, 1);
den   = 1 + (z.^2)./n_eff;
center= (p + (z.^2)./(2*n_eff)) ./ den;
half  = z .* sqrt( (p.*(1-p) + (z.^2)./(4*n_eff)) ./ n_eff ) ./ den;
lo = max(center - half, 0);
hi = min(center + half, 1);
lo(n==0) = NaN; hi(n==0) = NaN;
end

function s = g3(x)
if ~isfinite(x), s = 'NaN'; return; end
s = sprintf('%.3g', x);
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

function v = getfield2(M, f, default), if isempty(M) || ~isfield(M,f), v=default; else, v=M.(f); end, end
