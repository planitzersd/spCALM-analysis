function T = export_prism_xy_stagger(samples, y_fields, out_path, opts)
% EXPORT_PRISM_XY_STAGGER
% Build a Prism-friendly XY table with a single X column and "staggered"
% Y columns (NaN where a row does not belong to that sample/field).
%
% Inputs
%   samples  : struct array. For each s:
%              - samples(s).t        : Nx1 time vector (ordered)
%              - samples(s).<field>  : Nx1 fluorescence vector(s) (same length as t)
%              - samples(s).name     : (optional) string for labeling
%   y_fields : char or cellstr of fluorescence field names to export
%              e.g. 'x' or {'x','GFP','RFP'}
%   out_path : (optional) filename to write ('.csv', '.xlsx', or '.xls'); if empty, no file saved
%   opts     : (optional) struct:
%              .time_field   (default 't')
%              .add_nan_gap  (default true)   % insert a NaN separator row between samples
%              .gap_rows     (default 1)      % number of NaN rows between samples
%              .name_field   (default 'name') % where to read sample name
%              .prefix       (default '')     % prefix for Y column names
%              .overwrite    (default true)   % if Excel file exists, delete before writing
%              .sheet_base   (default 'Block')% base name for sheet chunks
%
% Output
%   T : table with variables:
%       - X : concatenated time
%       - one Y column per (sample × field), named "<sample>__<field>"

    if nargin < 2 || isempty(y_fields)
        error('Please provide y_fields (e.g., ''x'' or {''x'',''GFP''}).');
    end
    if ischar(y_fields) || isstring(y_fields)
        y_fields = cellstr(y_fields);
    end
    if nargin < 3
        out_path = '';
    end
    if nargin < 4
        opts = struct;
    end

    % defaults
    def.time_field  = 't';
    def.add_nan_gap = true;
    def.gap_rows    = 1;
    def.name_field  = 'name';
    def.prefix      = '';
    def.overwrite   = true;
    def.sheet_base  = 'Block';
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
    end

    S = numel(samples);
    assert(S >= 1, 'samples is empty.');

    % Helper to get a safe sample name
    function nm = sample_name(sidx)
        nm = sprintf('sample_%02d', sidx);
        if isfield(samples(sidx), opts.name_field) && ~isempty(samples(sidx).(opts.name_field))
            nm = string(samples(sidx).(opts.name_field));
        end
        nm = matlab.lang.makeValidName(char(nm));
    end

    % Pre-build column names for all (sample × field) Y columns
    Ynames = strings(S*numel(y_fields), 1);
    k = 1;
    for s = 1:S
        base = sample_name(s);
        for f = 1:numel(y_fields)
            Ynames(k) = matlab.lang.makeValidName(sprintf('%s%s__%s', opts.prefix, base, char(y_fields{f})));
            k = k + 1;
        end
    end

    % We'll build X as we go; for Y columns, we keep numeric arrays and expand as needed.
    X = [];                                  % concatenated time
    Ycols = cell(numel(Ynames), 1);          % cell of numeric vectors, each grows as X grows
    for i = 1:numel(Ycols), Ycols{i} = []; end

    % Iterate samples, append segment, and fill Y appropriately
    for s = 1:S
        % ---- fetch time ----
        tf = opts.time_field;
        assert(isfield(samples(s), tf), 'Sample %d missing time field "%s".', s, tf);
        t = samples(s).(tf);
        t = t(:);
        Ns = numel(t);
        if Ns == 0, continue; end

        % ---- build this segment’s X (plus optional NaN gap after) ----
        segX  = t;
        gapN  = (opts.add_nan_gap) * opts.gap_rows;
        if gapN > 0
            segX  = [segX; NaN(gapN,1)];
        end

        % ---- collect y vectors for this sample ----
        segY = cell(numel(y_fields), 1);
        for f = 1:numel(y_fields)
            fld = y_fields{f};
            assert(isfield(samples(s), fld), 'Sample %d missing field "%s".', s, fld);
            y = samples(s).(fld);
            y = y(:);
            assert(numel(y) == Ns, 'Length mismatch in sample %d: "%s" must match %s.', s, fld, tf);
            if gapN > 0
                y = [y; NaN(gapN,1)]; %#ok<AGROW>
            end
            segY{f} = y;
        end

        % ---- extend master X and all Y columns ----
        X = [X; segX];                 %#ok<AGROW>
        addL = numel(segX);

        baseIdx = (s-1)*numel(y_fields);
        for g = 1:numel(Ycols)
            if g > baseIdx && g <= baseIdx + numel(y_fields)
                fidx = g - baseIdx;
                Ycols{g} = [Ycols{g}; segY{fidx}]; %#ok<AGROW>
            else
                Ycols{g} = [Ycols{g}; NaN(addL,1)]; %#ok<AGROW>
            end
        end
    end

    % ---- assemble table ----
    T = table(X);
    for i = 1:numel(Ycols)
        T.(char(Ynames(i))) = Ycols{i};
    end

    % ---- write file if requested ----
    if ~isempty(out_path)
        [folder,~,ext] = fileparts(out_path);
        if ~isempty(folder) && ~exist(folder,'dir')
            mkdir(folder);
        end

        switch lower(ext)
            case '.csv'
                % CSV has no Excel row limits
                writetable(T, out_path);

            case {'.xlsx','.xlsm','.xls'}
                % Excel: enforce row limits and chunk across sheets
                switch lower(ext)
                    case {'.xlsx','.xlsm'}
                        maxRows = 1048576;
                    case '.xls'
                        maxRows = 65536;
                end

                if opts.overwrite && exist(out_path,'file')
                    delete(out_path);
                end

                N = height(T);
                if N <= maxRows
                    writetable(T, out_path, 'FileType','spreadsheet', 'Sheet', sprintf('%s01', opts.sheet_base));
                else
                    nBlocks = ceil(N / maxRows);
                    for b = 1:nBlocks
                        i1 = (b-1)*maxRows + 1;
                        i2 = min(b*maxRows, N);
                        Tb = T(i1:i2, :);
                        sheetName = sprintf('%s%02d', opts.sheet_base, b);
                        writetable(Tb, out_path, 'FileType','spreadsheet', 'Sheet', sheetName);
                    end
                end

            otherwise
                % default to CSV if unknown extension
                writetable(T, [out_path '.csv']);
        end
    end
end
