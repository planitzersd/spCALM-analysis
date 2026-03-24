function plot_select_samples(OUT, opts)
% PLOT_SELECT_SAMPLES
% Interactive viewer to overlay multiple sample trajectories with options:
% - Paged checkbox list (with slider) to pick which samples to plot
% - Choose model curve: 'logistic4', 'gamma', or 'none'
% - Toggle binned data error bars (binomial SE) or points-only
%
% INPUTS
%   OUT  : struct from fit_binned_kinetics* with OUT.samples(s) having:
%          .bins (table with t,p,n,t_lo,t_hi), .logistic4, .gamma, .name
%   opts : (optional) struct:
%          .visible_count   (default 16) #rows visible in checkbox list
%          .default_model   (default 'logistic4') 'logistic4'|'gamma'|'none'
%          .show_errorbars  (default true)
%          .show_points     (default true)
%          .line_width      (default 1.5) model curve width
%          .marker_size     (default 18) scatter marker size
%
% Example:
%   plot_select_samples(OUT);
%   plot_select_samples(OUT, struct('visible_count',20,'default_model','gamma'));

    if nargin < 2, opts = struct; end
    opts = set_default(opts, 'visible_count', 16);
    opts = set_default(opts, 'default_model', 'logistic4'); % 'gamma'|'none'
    opts = set_default(opts, 'show_errorbars', true);
    opts = set_default(opts, 'show_points', true);
    opts = set_default(opts, 'line_width', 1.5);
    opts = set_default(opts, 'marker_size', 18);

    S = OUT.samples;
    Ns = numel(S);
    if Ns == 0, error('OUT.samples is empty.'); end
    VISIBLE_COUNT = max(1, round(opts.visible_count));

    % Selection state
    selected = false(Ns,1);
    rowToIndex = nan(VISIBLE_COUNT,1);
    startIdx = 1;

    % ---------- UIFigure & top-level split ----------
    f = uifigure('Name','Kinetics Viewer','Position',[100 100 1180 720]);

    g = uigridlayout(f,[1 2]);
    g.ColumnWidth = {320, '1x'};
    g.RowHeight   = {'1x'};
    g.Padding     = [8 8 8 8];

    % Left panel (controls)
    leftPanel = uipanel(g, 'Title','Select Samples', 'FontWeight','bold');
    leftPanel.Layout.Row = 1;
    leftPanel.Layout.Column = 1;
    leftPanel.Scrollable = 'on';

    % Right panel (plot); disable autoresize so SizeChangedFcn runs
    rightPanel = uipanel(g, 'Title','Plot', 'FontWeight','bold');
    rightPanel.Layout.Row = 1;
    rightPanel.Layout.Column = 2;
    rightPanel.AutoResizeChildren = 'off';

    % Axes
    ax = uiaxes(rightPanel);
    ax.Box = 'on';
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Fraction fused (p)';
    ax.YLim = [0 1];
    grid(ax,'on');

    % Reposition axes when panel resizes
    rightPanel.SizeChangedFcn = @(~,~) positionAxes();

    % Initial position of axes
    positionAxes();

    % ---------- Left panel grid ----------
    % Rows: [header (70px), VISIBLE_COUNT rows of 22px, page label (fit), slider (30px)]
    lp = uigridlayout(leftPanel,[VISIBLE_COUNT+3 2]);
    lp.Padding     = [8 8 8 8];
    lp.RowHeight   = [ {70}, num2cell(repmat(22,1,VISIBLE_COUNT)), {'fit'}, {30} ];
    lp.ColumnWidth = {24, '1x'};

    % ----- Header controls (row 1, span 2 cols) -----
    hdr = uigridlayout(lp, [3 2]);
    hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];
    hdr.RowHeight   = {24, 24, 24};
    hdr.ColumnWidth = {'1x', '1x'};
    hdr.Padding = [0 0 0 0];

    % Model dropdown
    lblModel = uilabel(hdr,'Text','Model:','HorizontalAlignment','right');
    lblModel.Layout.Row = 1; lblModel.Layout.Column = 1;

    ddModel = uidropdown(hdr, 'Items', {'logistic4','gamma','none'}, ...
        'Value', lower(opts.default_model), ...
        'ValueChangedFcn', @(~,~) redrawPlot());
    ddModel.Layout.Row = 1; ddModel.Layout.Column = 2;

    % Error bars
    cbErr = uicheckbox(hdr,'Text','Error bars','Value', logical(opts.show_errorbars), ...
        'ValueChangedFcn', @(~,~) redrawPlot());
    cbErr.Layout.Row = 2; cbErr.Layout.Column = 1;

    % Points
    cbPts = uicheckbox(hdr,'Text','Points','Value', logical(opts.show_points), ...
        'ValueChangedFcn', @(~,~) redrawPlot());
    cbPts.Layout.Row = 2; cbPts.Layout.Column = 2;

    % Select all/none
    btnAll  = uibutton(hdr,'Text','Select All','ButtonPushedFcn', @(~,~) selectAll(true));
    btnNone = uibutton(hdr,'Text','Select None','ButtonPushedFcn', @(~,~) selectAll(false));
    btnAll.Layout.Row = 3;  btnAll.Layout.Column  = 1;
    btnNone.Layout.Row = 3; btnNone.Layout.Column = 2;

    % ----- Checkbox list (rows 2..VISIBLE_COUNT+1) -----
    chk = gobjects(VISIBLE_COUNT,1);
    lbl = gobjects(VISIBLE_COUNT,1);
    for r = 1:VISIBLE_COUNT
        chk(r) = uicheckbox(lp, 'Text','', 'Value',false, ...
            'ValueChangedFcn', @(h,~) onCheckChanged(h,r));
        chk(r).Layout.Row = r+1;  % header occupies row 1
        chk(r).Layout.Column = 1;

        lbl(r) = uilabel(lp, 'Text','', 'HorizontalAlignment','left');
        lbl(r).Layout.Row = r+1;
        lbl(r).Layout.Column = 2;
        lbl(r).WordWrap = 'on';
    end

    % ----- Page label (row VISIBLE_COUNT+2) -----
    lblPage = uilabel(lp, 'Text', '', 'HorizontalAlignment','left');
    lblPage.Layout.Row = VISIBLE_COUNT+2;
    lblPage.Layout.Column = 2;

    % ----- Slider (row VISIBLE_COUNT+3 spans 2 cols) -----
    sliderPanel = uipanel(lp, 'BorderType','none');
    sliderPanel.Layout.Row = VISIBLE_COUNT+3;
    sliderPanel.Layout.Column = [1 2];

    sld = uislider(sliderPanel, ...
        'Orientation','horizontal', ...
        'Limits', [1, max(1, Ns - VISIBLE_COUNT + 1)], ...
        'Value', 1, ...
        'MajorTicks', [], 'MinorTicks', []);
    sld.ValueChangingFcn = @(h,e) onSlide(round(e.Value));
    sld.ValueChangedFcn  = @(h,~) onSlide(round(h.Value));

    % Fit slider on resize
    sliderPanel.SizeChangedFcn = @(~,~) ...
        set(sld,'Position',[8, max(2,sliderPanel.Position(4)/2), max(40,sliderPanel.Position(3)-16), 3]);

    % Show slider only if needed
    sliderPanel.Visible = matlab.lang.OnOffSwitchState(Ns > VISIBLE_COUNT);

    % Initialize list + plot
    refreshList(startIdx);
    redrawPlot();

    % =============== nested helpers ===============

    function positionAxes()
        % Use inner position for padding-aware sizing
        ip = rightPanel.InnerPosition; % [x y w h]
        pad = 20;
        ax.Position = [ip(1)-15*pad, ip(2)+pad, max(50, ip(3)-2*pad), max(50, ip(4)-2*pad)];
    end

    function refreshList(firstIdx)
        startIdx = max(1, min(firstIdx, max(1, Ns - VISIBLE_COUNT + 1)));
        ending = min(Ns, startIdx + VISIBLE_COUNT - 1);

        for r = 1:VISIBLE_COUNT
            idx = startIdx + (r-1);
            if idx <= Ns
                rowToIndex(r) = idx;
                nm = '(unnamed)';
                if isfield(S(idx),'name') && ~isempty(S(idx).name)
                    nm = S(idx).name;
                end
                lbl(r).Text = sprintf('%s', nm);
                chk(r).Value = selected(idx);
                chk(r).Enable = 'on';
                lbl(r).Enable = 'on';
                chk(r).Visible = 'on';
                lbl(r).Visible = 'on';
            else
                rowToIndex(r) = NaN;
                lbl(r).Text = '';
                chk(r).Value = false;
                chk(r).Enable = 'off';
                lbl(r).Enable = 'off';
                chk(r).Visible = 'off';
                lbl(r).Visible = 'off';
            end
        end
        lblPage.Text = sprintf('Showing %d–%d of %d', startIdx, ending, Ns);
        sld.Limits = [1, max(1, Ns - VISIBLE_COUNT + 1)];
        sld.Value  = startIdx;
    end

    function onSlide(newStart)
        refreshList(newStart);
    end

    function onCheckChanged(h, r)
        idx = rowToIndex(r);
        if ~isnan(idx)
            selected(idx) = logical(h.Value);
            redrawPlot();
        end
    end

    function selectAll(flag)
        selected(:) = logical(flag);
        refreshList(startIdx);
        redrawPlot();
    end

    function redrawPlot()
        cla(ax); hold(ax,'on');
        whichModel = lower(ddModel.Value);
        showErr = logical(cbErr.Value);
        showPts = logical(cbPts.Value);

        selIdx = find(selected(:)');
        if isempty(selIdx)
            title(ax, 'No samples selected');
            drawnow; return;
        end

        cols = lines(max(1,numel(selIdx)));
        legends = {};

        for j = 1:numel(selIdx)
            sidx = selIdx(j);
            T = S(sidx).bins;
            if isempty(T) || height(T)==0, continue; end

            t  = T.t(:);
            p  = T.p(:);
            n  = T.n(:);
            se = sqrt(max(p.*(1-p)./max(n,eps), 0));

            % Data points / errors
            if showPts
                scatter(ax, t, p, opts.marker_size, ...
                    'MarkerEdgeColor', cols(j,:), 'MarkerFaceColor','none', ...
                    'DisplayName', sprintf('%s data', shortName(S(sidx))));
                legends{end+1} = sprintf('%s data', shortName(S(sidx))); %#ok<AGROW>
            end
            if showErr
                er = errorbar(ax, t, p, se, 'LineStyle','none', 'Color', cols(j,:)*0.55, ...
                    'DisplayName', sprintf('%s SE', shortName(S(sidx))));
                er.CapSize = 0;
            end

            % Model curve
            switch whichModel
                case 'logistic4'
                    if hasfield_chain(S(sidx), {'logistic4','params'})
                        prm = S(sidx).logistic4.params;
                        tt = linspace(min(t), max(t), 400);
                        ph = logistic4_eval(tt, prm.B, prm.A, prm.t50, prm.k);
                        plot(ax, tt, ph, '--', 'Color', cols(j,:), 'LineWidth', opts.line_width, ...
                            'DisplayName', sprintf('%s (4PL)', shortName(S(sidx))));
                        legends{end+1} = sprintf('%s 4PL', shortName(S(sidx))); %#ok<AGROW>
                    end
                case 'gamma'
                    if hasfield_chain(S(sidx), {'gamma','params'})
                        prm = S(sidx).gamma.params;
                        tt = linspace(min(t), max(t), 400);
                        ph = gamma_onset_eval(tt, prm.B, prm.A, prm.t0, prm.N, prm.k);
                        plot(ax, tt, ph, '--', 'Color', cols(j,:), 'LineWidth', opts.line_width, ...
                            'DisplayName', sprintf('%s (Gamma)', shortName(S(sidx))));
                        legends{end+1} = sprintf('%s Gamma', shortName(S(sidx))); %#ok<AGROW>
                    end
                otherwise
                    % none
            end
        end

        xlabel(ax,'Time (s)'); ylabel(ax,'Fraction fused (p)');
        ylim(ax,[0 1]); grid(ax,'on');
        if ~isempty(legends)
            legend(ax,'Location','bestoutside');
        end
        title(ax, sprintf('Selected: %d / %d', numel(selIdx), Ns));
        hold(ax,'off');
        drawnow;
    end

    function s = shortName(Ss)
        if isfield(Ss,'name') && ~isempty(Ss.name), s = char(Ss.name);
        else, s = 'sample'; end
    end

    % ---- local evaluators (copy from your fit code) ----
    function p = logistic4_eval(t, B, A, t50, k)
        p = B + A ./ (1 + exp(-k .* (t - t50)));
        p = min(max(p, 1e-9), 1 - 1e-9);
    end
    function p = gamma_onset_eval(t,B,A,t0,N,k)
        tau = max(t - t0, 0);
        theta = 1 / max(k, eps);
        p = B + A .* gamcdf(tau, max(N,1e-6), theta);
        p = min(max(p, 1e-9), 1 - 1e-9);
    end

    % ---- small utils ----
    function tf = hasfield_chain(s, chain)
        tf = true;
        for ii = 1:numel(chain)
            f = chain{ii};
            if ~isfield(s, f) || isempty(s.(f)), tf = false; return; end
            s = s.(f);
        end
    end
    function s = set_default(s, field, val)
        if ~isfield(s, field) || isempty(s.(field)), s.(field) = val; end
    end
end
