%========================================================================%
%~~~~~~~~~~~~~~ BEGIN MAIN, SPECIFY USER DEFINED PARAMETERS ~~~~~~~~~~~~~%
%========================================================================%

%USER OPTIONS%

% where to save output
outputfoldername = 'chirpfit11';

% folder to open to to look for .csv data files
homepath = '';

%Annotation tag order ('_'-separated, excluding first clause)
tag_fields = {'typetag','recsattag', 'NAItag', 'pHtag', 'reptag'};
sort_order = [5, 4, 3, 2, 1]; %sort tags in this order
istagnumeric = [false, false, false, false, false]; %treat tag as number?

% Choose plateau guesses based on flowjo data
base_sec = 10; %time (from run start) to consider for pre-pH-drop averaging
plat_sec = 60; %time (from run end) to consider for post-pH-drop averaging

%Use a dialog to choose samples to base unfused and fused states on, these
%should be samples that reach near maximal fusion. Or, specify their sample
%indices. Or, don't specify anything and compute on all samples.

opts.emissions_use_dialog = true;
%opts.emissions_selected_samples = [23,24,25,33,34,35];

%Prespecify ranges of time in seconds for computing unfused/fused
%probabilities based on selected samples above, or just leave empty for
%auto selection. 
opts.emissions_plateau_cutoffs = [base_sec plat_sec];
%opts.emissions_plateau_cutoffs = [];

% You may need to tune other parameters in the adaptive bin, compute frac 
% fused, fit models section.
%%
%=========================================================================%
%~~~~~~~~~~~~~~~~~~ IMPORT CSV DATA (FLOWJO OUTPUT) ~~~~~~~~~~~~~~~~~~~~~~%
%=========================================================================%

%Select .csv files for analysis
[files1, paths1] = uigetfile('*.csv', 'Select .csv files exported from FlowJo', homepath, "MultiSelect", "on");
target_folder = fullfile(paths1, outputfoldername);
%%
%If there is a corresponding free virus .csv, write the suffix of the
%filename preceding .csv below
file1suffix = ['Pairs'];
%file2suffix = 'FreeVirus';

%files2 = replace(files1,file1suffix,file2suffix); %get filenames of corresponding .csvs for other gates

if iscell(files1)
        kkmax = length(files1);
else
        kkmax = 1;
end

for kk=1:kkmax
    %try
    %%fraqDeq = [];
    nbinevents = [];
    medRFPA = [];
    medFITCA = [];
    medAPCA = [];
    medVSSCH = [];

    if iscell(files1)
        filestemp1 = files1{kk};
        %filestemp2 = files2{kk};
    else
        filestemp1 = files1;
        %filestemp2 = files2;
    end
    
    sprintf('Loading %d/%d: %s', kk,kkmax,filestemp1)
    [Cevents, runduration,filename] = loadgatedcsv(paths1,filestemp1);
    %[Ceventsfreevirus, ~, ~] = loadgatedcsv(paths1,filestemp2);
   
    
    output(kk).filename = filename;
    output(kk).path = paths1;

    RFP = str2double(Cevents.RFP_A(1:end-1));
    FITC = str2double(Cevents.FITC_H(1:end-1));
    VSSC = str2double(Cevents.VioletSSC_H(1:end-1));
    APC = str2double(Cevents.APC_A(1:end-1));
    time = str2double(Cevents.Time(1:end-1));
    time = time.*runduration./max(time,[],'all');

    totaltimemin = max(time);

    samples(kk).t = time;

    %samples(kk).x = ltObj.transform(abs(RFP));
    samples(kk).x = asinh(RFP);
    samples(kk).name = filename;
    
    rawsamples(kk).t = time;
    rawsamples(kk).x = RFP;

    samples(kk).RFP = RFP;
    samples(kk).FITC = FITC;
    samples(kk).APC = APC;
    samples(kk).FITCtoAPCratio = FITC./APC;
    samples(kk).VSSC = VSSC;
end

%%
%=========================================================================%
%~~~~~~~~~~~~ ADAPTIVE BIN, COMPUTE FRAC FUSED, FIT MODELS ~~~~~~~~~~~~~~~%
%=========================================================================%
opts.maxbins = 5000;
opts.target_count_mode = 'by-total';
opts.target_count_scale = 0.001;   % e.g., 1% of events per bin
opts.target_count_min   = 50;     % optional tighter clamps
opts.target_count_max   = 1000;

opts.display = 'off';
opts.computeCI = true;
opts.pre_t0_bins = 2;

opts.logistic_weight_mode   = 'slope';   % or 'slope' or ''
opts.logistic_optimize_beta = false;      % turn on R^2-based beta selection
%opts.logistic_beta_grid     = linspace(0,10,11); % can expand/tighten
%opts.logistic_beta_grid     = [5 6 7 8 9 10]; % can expand/tighten
opts.logistic_weight_beta = 0;

opts.t0_from_rate = true;   % to populate OUT.samples(s).t0
opts.gamma_t0_mode = 'prior';
opts.gamma_t0_prior_scale = 0.5;

opts.gamma_weight_mode   = 'slope-logistic';   % or 'slope' or ''
opts.gamma_optimize_beta = false;      % turn on R^2-based beta selection
%opts.gamma_beta_grid     = linspace(0,10,11); % can expand/tighten
%opts.gamma_beta_grid     = [5 6 7 8 9 10]; % can expand/tighten
opts.gamma_weight_beta = 0;


opts.weight_floor = 20; 

opts.gamma_init_quants = [0.2 0.4 0.6 0.8];
opts.gamma_init_weights = ones(size(opts.gamma_init_quants));
%opts.gamma_init_weights = [1.5 1.5 1.5 1.5];

OUT = fit_binned_kinetics12(samples, opts);

%plot_binned_kinetics(OUT, 10)
% Learned global emissions:
% OUT.mu1, OUT.mu2, OUT.sigma

% % Per-sample final fits with SEs and 95% CIs:
% OUT.samples(s).final.params
% OUT.samples(s).final.se
% OUT.samples(s).final.ci
% OUT.samples(s).final.bic

% indx = 55;

% Example: visualize logistic weighting for sample 2
% plot_transition_weights(OUT.samples(indx).logistic4, ...
%     'Time (s)', 'Proportion p', OUT.samples(indx).name);

%Example: visualize gamma fit weights
% plot_transition_weights(OUT.samples(indx).gamma, ...
%     'Time (s)', 'Proportion p', OUT.samples(indx).name);

%%

OUT.samples = annotatestructure(OUT.samples, tag_fields, sort_order, istagnumeric);
samples = annotatestructure(samples, tag_fields, sort_order, istagnumeric);

%%
%=========================================================================%
%~~~~~~~~~~~~~~~~~~~~~ PROCESS FITTED MODEL RESULTS ~~~~~~~~~~~~~~~~~~~~~~%
%=========================================================================%

%organize fit outputs into tables for pasting into Graphpad Prism

tags = {tag_fields{1:end-1}};

model_name   = 'logistic4';                          % or 'logistic4' / 'gamma'
param_fields = {'span', 't50', 'k'};                        % columns per replicate

[T4PLparam,T4PLgof] = make_prism_table(OUT, tag_fields, model_name, param_fields, target_folder);

model_name   = 'gamma';                          % or 'logistic4' / 'gamma'
param_fields = {'span', 'N', 'k', 't0'};   
[Tgammaparam,Tgammagof] = make_prism_table(OUT, tag_fields, model_name, param_fields, target_folder);

opts = struct('t0_se_mode','auto', ...  % or 'fixed' with t0_se_value, or 'zero'
              'alpha', 0.05);

[T_samp, T_cond] = summarize_delta_t50_t0(OUT, tag_fields, target_folder, opts);
[T_deviationoft0ingamfit_samp,T_deviationoft0ingamfit_comp] = summarize_t0_shift(OUT, tag_fields, target_folder, ...
    struct('alpha',0.05, 'den_floor',1, 'ratio_mode','rate', ...
           'excel_name','t0_shift.xlsx', 'csv_prefix','t0_shift'));

[Tlog_samp, Tlog_cond, Tgam_samp, Tgam_cond, Tbest_samp, Tbest_cond] = ...
    summarize_fit_parameters(OUT, tag_fields, target_folder);

%%
%Perform very simple binning of fluorescence parameters

[~, edges] = histcounts(samples(1).t,300);
centers = (edges(2:end)+edges(1:end-1))/2;
FITCtoAPCratioBINNED = nan(length(centers),length(samples));
RFPBINNED = nan(length(centers),length(samples));
COUNTSBINNED = nan(length(centers),length(samples));

for kk=1:length(samples)
    curtime = samples(kk).t;
    currat = samples(kk).FITCtoAPCratio;
    currfp = samples(kk).RFP;
    for jj = 1:length(centers)
        indices = find(curtime >= edges(jj) & curtime <= edges(jj+1));
        FITCtoAPCratioBINNED(jj,kk) = median(currat(indices));
        RFPBINNED(jj,kk) = median(currfp(indices));
        COUNTSBINNED(jj,kk) = length(indices);
    end
end

%%
% Plot adaptive-binned trajectories with model fits, save to pdf
pdf_file = fullfile(target_folder, ['kinetics_report' num2str(ceil(1000*rand)) '.pdf']);

plot_binned_kinetics_pdf(OUT, pdf_file, struct('per_page', 2));
%%
% Interactive plotting tool: compare specified trajectories
plot_select_samples(OUT)

%%
%=========================================================================%
%~~~~~~~~~~~~~~~~COMPUTE MODEL-AGNOSTIC SPAN AND LAGTIME ~~~~~~~~~~~~~~~~~%
%=========================================================================%

opts = struct('base_sec', base_sec, 'plat_sec', plat_sec);

% Compute span
opts_span = opts;
opts_span.excel_path = fullfile(target_folder, 'span_simple.xlsx');
T_span = span_simple(OUT, opts_span);



% Compute t50
opts_tq = opts;
opts_tq.q = 0.5;

for kk = 1:length(OUT.samples)
    tempstruct(kk).t = OUT.samples(kk).bins.t;
    tempstruct(kk).p = OUT.samples(kk).bins.p;
    tempstruct(kk).name = OUT.samples(kk).name;
end

for kk = 1:length(OUT.samples)
    tempstruct(kk).t = OUT.samples(kk).bins.t;
    tempstruct(kk).r = OUT.samples(kk).bins.p;
    tempstruct(kk).name = OUT.samples(kk).name;
    tempstruct(kk).t_lo = OUT.samples(kk).bins.t_lo;
    tempstruct(kk).t_hi = OUT.samples(kk).bins.t_hi;
end

%Acquire median effect times
alt_t50_est = tq_isotonic_struct_2(tempstruct, opts_tq);
alt_t50_est = annotatestructure(alt_t50_est, tag_fields, sort_order, istagnumeric);

%%
% Plot median effect times and spans over individual trajectories 

plot_opts = struct('show_span', true, 'res_span', T_span, ...
                   'show_windows', true, ...
                   'base_sec', base_sec, 'plat_sec', plat_sec, ...
                   'save_path', fullfile(target_folder, '/span_lagtime_plots/'), ...
                   'close_after_save', true, ...
                   'save_format', 'png');

plot_tq_check(OUT, alt_t50_est, [], plot_opts);

%%

%compute lag times (tq-t0)
conditiontags = tag_fields(1:3);

ivw_t50_est = aggregate_tq_by_condition_ivw(alt_t50_est, conditiontags, 0.05);

[ivw_lagtime_cond, ivw_lagtime_rep] = subtract_baseline_and_pool_ivw(alt_t50_est, OUT.samples, conditiontags, 0.05);

%save lagtime tables
writetable(ivw_lagtime_cond, fullfile(target_folder, 'ivw_lagtime_cond.xlsx'));
writetable(ivw_lagtime_rep, fullfile(target_folder, 'ivw_lagtime_rep.xlsx'));


%%
% organize adaptive-binned frac. fused trajectories for plotting in prism
staggeredBinData = export_prism_xy_stagger(tempstruct, {'p'});

%%

%save whole workspace for later work
save(fullfile(target_folder, 'workspace.mat'));


%========================================================================%
%~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%========================================================================%

%%

%========================================================================%
%~~~~~~~~~~~~~~~~~~~~~~~~ HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%========================================================================%

function [eventtable, runduration, filename] = loadgatedcsv(path, file)

    path = [path file];

    pathsplit = split(path,'\');
    filename = split(pathsplit{end},'.csv');
    filename = filename{1};

    Ckeywords = readcell(path, "Range", "A1:B1000");
    
    btimidx = find(strcmp(Ckeywords,'$BTIM'));
    etimidx = find(strcmp(Ckeywords,'$ETIM'));
    btim = seconds(Ckeywords{btimidx,2});
    etim = seconds(Ckeywords{etimidx,2});
    runduration = etim-btim;
        
    varlineidx = find(strcmp(Ckeywords,'FSC-A'));

    opts = delimitedTextImportOptions('VariableNamesLine' , varlineidx(1) , ...
         'DataLines' , [varlineidx(1)+1 Inf] );
    
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames'); 
    eventtable = readtable(path,opts); 

end

function out = annotatestructure(in, tags, sortOrder, isNumeric)
    % ANNOTATESTRUCTURE Adds tag fields to structure array and sorts
    %
    % Inputs:
    %   in        - Structure array with 'name' field (format: "export_tag1_tag2_..." or "tag1_tag2_...")
    %   tags      - Cell array of tag names in order they appear in name
    %   sortOrder - Vector of indices specifying sort priority (e.g., [4,3,2,1])
    %   isNumeric - (Optional) Logical vector where 1=numeric, 0=char (default: all char)
    %
    % Output:
    %   out       - Annotated and sorted structure array
    
    % Set default: all tags are char unless specified
    if nargin < 4 || isempty(isNumeric)
        isNumeric = false(1, length(tags));
    end
    
    for jj = 1:length(in)
        try
           curname = in(jj).name;
           splname = strsplit(curname, '_');
            % Check if first substring is "export" and skip it
            startIdx = 1;
            if ~isempty(splname) && strcmpi(splname{1}, 'export')
                startIdx = 2;
            end
            
            % Assign each tag value to corresponding field
            for kk = 1:length(tags)
                tagIdx = startIdx + kk - 1;
                if tagIdx <= length(splname)
                    if isNumeric(kk)
                        % Convert to numeric
                        in(jj).(tags{kk}) = str2double(splname{tagIdx});
                    else
                        % Keep as char
                        in(jj).(tags{kk}) = char(splname{tagIdx});
                    end
                else
                    % Not enough substrings - set to empty
                    if isNumeric(kk)
                        in(jj).(tags{kk}) = NaN;
                    else
                        in(jj).(tags{kk}) = char();
                    end
                end
            end
            
        catch
            % If any error occurs, set all tags to empty char or NaN
            for kk = 1:length(tags)
                if isNumeric(kk)
                    in(jj).(tags{kk}) = NaN;
                else
                    in(jj).(tags{kk}) = char();
                end
            end
        end
    end
    
    % Convert structure array to table
    T = struct2table(in);
    
    % Sort by tags in the order specified by sortOrder
    % sortOrder is in reverse priority (last element is highest priority)
    for ii = 1:length(sortOrder)
        tagIdx = sortOrder(ii);
        if tagIdx >= 1 && tagIdx <= length(tags)
            T = sortrows(T, tags{tagIdx});
        end
    end
    
    % Convert back to structure array
    out = table2struct(T);
end


function out = annotatestructure_old(in)
    for jj = 1:length(in)
        curname = in(jj).name;
    
        splname = strsplit(curname,'_');

        try
            %if str2double(splname{2})==3 || str2double(splname{2})==6 || str2double(splname{2})==0
                in(jj).muttag = char(splname{1});
                in(jj).abtag = char(splname{2});
                in(jj).NAItag = char(splname{3});
                in(jj).reptag = str2double(splname{4});
                
                % in(jj).rectypetag = char(splname{2});
                % in(jj).recsattag = str2double(splname{3});
                % in(jj).NAItag = char(splname{4});
                % in(jj).reptag = char(splname{5});
        catch
                in(jj).muttag = char();
                in(jj).abtag = char();
                in(jj).NAItag = char();
                in(jj).reptag = char();

                % in(jj).rectypetag = char();
                % in(jj).recsattag = str2double('');
                % in(jj).NAItag = char();
                % in(jj).reptag = char();
        end
    end
    
    
    % Convert the structure array to a table
    T = struct2table(in);
    
    % Sort the table by the desired field (e.g., 'score')

    % sortedT = sortrows(T, 'reptag'); % Ascending order
    % sortedT = sortrows(sortedT, 'pHtag'); % Ascending order
    % sortedT = sortrows(sortedT, 'NAItag'); % Ascending order
    % sortedT = sortrows(sortedT, 'recsattag'); % Ascending order
    
    sortedT = sortrows(T, 'reptag'); % Ascending order
    sortedT = sortrows(sortedT, 'abtag'); % Ascending order
    sortedT = sortrows(sortedT, 'NAItag'); % Ascending order
    sortedT = sortrows(sortedT, 'muttag'); % Ascending order
    
    % For descending order: sortedT = sortrows(T, 'score', 'descend');
    
    % Convert the sorted table back to a structure array (optional)
    sortedS = table2struct(sortedT);
    out = sortedS;
end



function [Tparam, Tgof] = make_prism_table(OUT, tag_fields, model_name, param_fields, target_folder, effect_q)
% MAKE_PRISM_TABLE
% Assemble parameter and goodness-of-fit tables in Prism-friendly wide format,
% grouped by user-provided tag fields and split across replicate columns.
%
% Adds (for gamma model): per-replicate columns showing:
%   - t_reach_<q>_gamma__repX  : time to reach effect_q fraction (default 0.99)
%   - reached_<q>_gamma__repX  : boolean (1/0) whether t_reach <= run duration
%
% Inputs:
%   OUT          : struct with OUT.samples(s) having .<model_name> fits and .bins
%   tag_fields   : cellstr of tag field names present on OUT.samples(s)
%                  (e.g., {'pHtag','NAItag','recsattag','reptag'})
%   model_name   : 'logistic4' or 'gamma'
%   param_fields : cellstr of model parameter names to export
%                  (e.g., {'span','bottom','t50','k'} for 4PL; {'span','bottom','t0','N','k'} for gamma)
%   target_folder: folder to save CSVs (optional; '' or [] to skip saving)
%   effect_q     : (optional) fraction for gamma reach columns; default 0.99
%
% Outputs:
%   Tparam : wide parameter table (one row per condition; columns per replicate)
%   Tgof   : wide GOF table (R2/AIC/BIC per replicate)

if nargin < 6 || isempty(effect_q), effect_q = 0.99; end
if nargin < 5, target_folder = ''; end
if ischar(tag_fields) || isstring(tag_fields), tag_fields = cellstr(tag_fields); end
if ischar(param_fields) || isstring(param_fields), param_fields = cellstr(param_fields); end

S = OUT.samples;
Ns = numel(S);
assert(Ns>=1, 'No samples in OUT.samples');

% --------- Normalize tag fields to strings for grouping ----------
tags_tbl = table();
for i = 1:numel(tag_fields)
    f = tag_fields{i};
    tags_tbl.(f) = repmat(string(missing), Ns, 1);
    for s = 1:Ns
        if isfield(S(s), f) && ~isempty(S(s).(f))
            val = S(s).(f);
            if isstring(val) || ischar(val)
                tags_tbl.(f)(s) = string(val);
            elseif isnumeric(val) || islogical(val)
                tags_tbl.(f)(s) = string(val);
            else
                try
                    tags_tbl.(f)(s) = string(val);
                catch
                    tags_tbl.(f)(s) = "<unk>";
                end
            end
        end
    end
end

% Identify replicate field if present
has_reptag = any(strcmpi(tag_fields, 'reptag'));

% Build condition key (excluding reptag)
key_fields = tag_fields(~strcmpi(tag_fields, 'reptag'));
if isempty(key_fields)
    cond_key = repmat("ALL", Ns, 1);
else
    cond_key = join(tags_tbl{:, key_fields}, "__", 2);
end
% normalize missing keys
cond_key = string(cond_key);
cond_key(ismissing(cond_key)) = "MISSING";
[uniq_keys, ~, key_idx] = unique(cond_key);

% --------- Collect per-sample model params & GOF & run duration ----------
MNAME = model_name; % alias
SampleParam = cell(Ns, 1);
SampleGOF   = cell(Ns, 1);
SampleRunEnd = nan(Ns,1);

for s = 1:Ns
    % Run duration (from bins)
    run_end = NaN;
    if isfield(S(s),'bins') && ~isempty(S(s).bins)
        T = S(s).bins;
        if ismember('t_hi', T.Properties.VariableNames)
            run_end = max(T.t_hi);
        elseif ismember('t', T.Properties.VariableNames)
            run_end = max(T.t);
        end
    end
    SampleRunEnd(s) = run_end;

    % Extract params
    P = struct();
    if isfield(S(s), MNAME) && ~isempty(S(s).(MNAME)) && isfield(S(s).(MNAME),'params')
        prm = S(s).(MNAME).params;
        for j = 1:numel(param_fields)
            pf = param_fields{j};
            P.(pf) = getfield_soft(prm, pf, NaN);
        end
        % convenience extras
        if ~isfield(P,'span'),   P.span   = getfield_soft(prm, 'span', NaN); end
        if ~isfield(P,'bottom'), P.bottom = getfield_soft(prm, 'bottom', NaN); end
        if ~isfield(P,'top'),    P.top    = getfield_soft(prm, 'top', NaN); end
    else
        for j = 1:numel(param_fields), P.(param_fields{j}) = NaN; end
        P.span = NaN; P.bottom = NaN; P.top = NaN;
    end
    SampleParam{s} = P;

    % Compute GoF
    SampleGOF{s} = compute_gof(S(s), MNAME);
end

% --------- If gamma model: precompute t_reach and "reached" flag ----------
tReach = nan(Ns,1);
reached = nan(Ns,1);
if strcmpi(MNAME, 'gamma')
    for s = 1:Ns
        P = SampleParam{s};
        N = getfield_soft(P, 'N', NaN);
        k = getfield_soft(P, 'k', NaN);
        t0= getfield_soft(P, 't0', NaN);
        if all(isfinite([N,k,t0])) && N > 0 && k > 0
            scale = 1/max(k,eps); % MATLAB's gaminv uses shape, scale
            try
                tReach(s) = gaminv(effect_q, N, scale) + t0;
            catch
                tReach(s) = NaN;
            end
            if isfinite(SampleRunEnd(s))
                reached(s) = double(tReach(s) <= SampleRunEnd(s));
            else
                reached(s) = NaN;
            end
        else
            tReach(s) = NaN; reached(s) = NaN;
        end
    end
end

% --------- Build replicate labels per sample ----------
if has_reptag
    rep_vals = tags_tbl.('reptag');
else
    % Create pseudo replicate IDs within each condition (encounter order)
    rep_vals = strings(Ns,1);
    counterMap = containers.Map('KeyType','char','ValueType','double');
    for s = 1:Ns
        kstr = uniq_keys(key_idx(s));
        if ismissing(kstr), kstr = "MISSING"; end
        keyChar = char(kstr);
        if ~isKey(counterMap, keyChar), counterMap(keyChar) = 0; end
        counterMap(keyChar) = counterMap(keyChar) + 1;
        rep_vals(s) = string(counterMap(keyChar));
    end
end
rep_vals(ismissing(rep_vals)) = "MISSING";

% ---- GLOBAL replicate roster (so every row has the same columns) ----
if has_reptag
    global_reps = unique(rep_vals);
    global_reps(ismissing(global_reps)) = "MISSING";
else
    % no reptag: size to max #samples per condition
    R = 0;
    for k = 1:numel(uniq_keys)
        R = max(R, sum(key_idx == k));
    end
    global_reps = string(1:R);
end

% --------- Build wide tables per condition with GLOBAL replicate columns ----------
ParamRows = {};
GOFRows   = {};

for k = 1:numel(uniq_keys)
    key = uniq_keys(k);
    if ismissing(key), key = "MISSING"; end

    % Base tag columns (keys)
    baseRow = table();
    for i = 1:numel(key_fields)
        f = key_fields{i};
        rows = find(key_idx == k);
        val = tags_tbl.(f)(rows(1));
        if ismissing(val), val = "MISSING"; end
        baseRow.(f) = val;
    end

    rows_k = find(key_idx == k);  % samples in this condition

    % ---------- PARAM row ----------
    prow = baseRow;

    % Pre-create parameter columns for ALL global replicates
    for j = 1:numel(param_fields)
        pf = param_fields{j};
        for r = 1:numel(global_reps)
            repval = global_reps(r);
            colname = matlab.lang.makeValidName(sprintf('%s__rep%s', pf, char(repval)));
            prow.(colname) = NaN;
        end
    end

    % Gamma extras (time to reach q and reached flag) – only when exporting gamma
    add_gamma_extras = strcmpi(MNAME,'gamma');
    if add_gamma_extras
        tagq = strrep(sprintf('%.2g',effect_q),'.','p'); % e.g., 0.99 -> 0p99
        for r = 1:numel(global_reps)
            repval = global_reps(r);
            col1 = matlab.lang.makeValidName(sprintf('t_reach_%s_gamma__rep%s', tagq, char(repval)));
            col2 = matlab.lang.makeValidName(sprintf('reached_%s_gamma__rep%s', tagq, char(repval)));
            prow.(col1) = NaN;
            prow.(col2) = NaN;
        end
    end

    % Fill the replicate slots for the samples in this condition
    for idxRow = 1:numel(rows_k)
        s = rows_k(idxRow);
        repval = rep_vals(s);
        if ismissing(repval), repval = "MISSING"; end
        ridx = find(global_reps == repval, 1, 'first');
        if isempty(ridx), continue; end

        % params
        P = SampleParam{s};
        for j = 1:numel(param_fields)
            pf = param_fields{j};
            val = getfield_soft(P, pf, NaN);
            colname = matlab.lang.makeValidName(sprintf('%s__rep%s', pf, char(repval)));
            prow.(colname) = val;
        end

        % gamma extras
        if add_gamma_extras
            tagq = strrep(sprintf('%.2g',effect_q),'.','p');
            col1 = matlab.lang.makeValidName(sprintf('t_reach_%s_gamma__rep%s', tagq, char(repval)));
            col2 = matlab.lang.makeValidName(sprintf('reached_%s_gamma__rep%s', tagq, char(repval)));
            prow.(col1) = tReach(s);
            prow.(col2) = reached(s);
        end
    end

    ParamRows{end+1} = prow; %#ok<AGROW>

    % ---------- GOF row ----------
    grow = baseRow;
    gof_fields = {'R2','AIC','BIC'};
    % Pre-create GOF columns for ALL global replicates
    for g = 1:numel(gof_fields)
        gf = gof_fields{g};
        for r = 1:numel(global_reps)
            repval = global_reps(r);
            colname = matlab.lang.makeValidName(sprintf('%s__rep%s', gf, char(repval)));
            grow.(colname) = NaN;
        end
    end

    % Fill per replicate
    for idxRow = 1:numel(rows_k)
        s = rows_k(idxRow);
        repval = rep_vals(s);
        if ismissing(repval), repval = "MISSING"; end
        ridx = find(global_reps == repval, 1, 'first');
        if isempty(ridx), continue; end
        G = SampleGOF{s};
        for g = 1:numel(gof_fields)
            gf = gof_fields{g};
            colname = matlab.lang.makeValidName(sprintf('%s__rep%s', gf, char(repval)));
            grow.(colname) = getfield_soft(G, gf, NaN);
        end
    end

    GOFRows{end+1} = grow; %#ok<AGROW>
end

% Concatenate rows
Tparam = vertcat(ParamRows{:});
Tgof   = vertcat(GOFRows{:});

% Sort rows by tag_fields (stable)
if ~isempty(key_fields)
    Tparam = sortrows(Tparam, key_fields);
    Tgof   = sortrows(Tgof, key_fields);
end

% Save CSVs if requested
if ~isempty(target_folder)
    if ~exist(target_folder, 'dir'), mkdir(target_folder); end
    writetable(Tparam, fullfile(target_folder, sprintf('prism_params_%s.csv', lower(MNAME))));
    writetable(Tgof,   fullfile(target_folder, sprintf('prism_gof_%s.csv',    lower(MNAME))));
    if strcmpi(MNAME,'gamma')
        tagq = strrep(sprintf('%.2g',effect_q),'.','p'); % e.g., 0.99 -> 0p99
        writetable(Tparam, fullfile(target_folder, sprintf('prism_params_gamma_q%s.csv', tagq)));
    end
end
end

function [T_samples, T_conditions] = summarize_t0_shift(OUT, tag_fields, target_folder, opts)
% SUMMARIZE_T0_SHIFT
% Compare event-rate t0 vs. gamma-fit t0 per sample and by condition.
%
% Relative difference is defined (by default) as:
%    rel = (t0_gamma - t0_rate) / max(|t0_rate|, den_floor)
% and also reported as percent (100*rel).
%
% Inputs
%   OUT           : struct from fit_binned_kinetics*
%                   expects OUT.samples(s).t0, .t0_se, and .gamma.params.t0
%                   and (when available) .gamma.se.t0, .gamma.notes.used_fixed_t0
%   tag_fields    : cellstr of tag field names (e.g., {'recsattag','NAItag','rectypetag','reptag'})
%   target_folder : '' to skip writing, or a folder path to save CSV/Excel
%   opts (struct) :
%       .alpha        : CI tail probability (default 0.05 => 95% CI)
%       .den_floor    : floor for the relative denominator (seconds) (default 1)
%       .ratio_mode   : 'rate' (default) — denominator is |t0_rate|
%                       'abs'  — denominator is |t0_rate| (same as 'rate'; alias)
%                       'none' — report only absolute difference (rel = NaN, percent = NaN)
%       .excel_name   : filename for Excel (default 't0_shift.xlsx') if target_folder provided
%       .csv_prefix   : prefix for CSV base name (default 't0_shift')
%
% Outputs
%   T_samples    : per-sample table (t0_rate, t0_gamma, SEs, diff, relative, CI, flags, tags)
%   T_conditions : per-condition inverse-variance weighted mean of the absolute/relative shift
%
% Notes on uncertainty:
%   - SE(diff) is propagated as sqrt(SE_rate^2 + SE_gamma^2) when both available.
%   - For the relative metric, we divide SE(diff) by the (floored) denominator
%     (i.e., we ignore denominator uncertainty; set .ratio_mode='none' if undesired).

if nargin < 4, opts = struct; end
if nargin < 3, target_folder = ''; end
if ischar(tag_fields) || isstring(tag_fields), tag_fields = cellstr(tag_fields); end

% ---- defaults
dd.alpha      = 0.05;
dd.den_floor  = 1;             % seconds
dd.ratio_mode = 'rate';        % 'rate'|'abs'|'none'
dd.excel_name = 't0_shift.xlsx';
dd.csv_prefix = 't0_shift';
fn = fieldnames(dd);
for i=1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = dd.(fn{i}); end
end
z = norminv(1 - opts.alpha/2);

S  = OUT.samples;
Ns = numel(S);

% ---- pull per-sample values
Name     = strings(Ns,1);
t0_rate  = nan(Ns,1);
se_rate  = nan(Ns,1);
t0_gam   = nan(Ns,1);
se_gam   = nan(Ns,1);
used_fixed_t0 = false(Ns,1);

for s=1:Ns
    % name
    if isfield(S(s),'name') && ~isempty(S(s).name)
        Name(s) = string(S(s).name);
    else
        Name(s) = "sample_"+s;
    end

    % event-rate t0
    if isfield(S(s),'t0') && isfinite(S(s).t0)
        t0_rate(s) = S(s).t0;
        % SE for rate t0
        if isfield(S(s),'t0_se') && isfinite(S(s).t0_se)
            se_rate(s) = max(0, S(s).t0_se);
        elseif isfield(OUT,'opts') && isfield(OUT.opts,'t0_binwidth') && isfinite(OUT.opts.t0_binwidth)
            se_rate(s) = OUT.opts.t0_binwidth / sqrt(12);
        else
            se_rate(s) = NaN;
        end
    end

    % gamma t0 (optimized or fixed)
    if isfield(S(s),'gamma') && ~isempty(S(s).gamma) && isfield(S(s).gamma,'params') ...
            && isfield(S(s).gamma.params,'t0') && isfinite(S(s).gamma.params.t0)
        t0_gam(s) = S(s).gamma.params.t0;

        if isfield(S(s).gamma,'notes') && isfield(S(s).gamma.notes,'used_fixed_t0')
            used_fixed_t0(s) = logical(S(s).gamma.notes.used_fixed_t0);
        end

        % SE for gamma t0: 0 if fixed; otherwise from fit if present
        if used_fixed_t0(s)
            se_gam(s) = 0;
        elseif isfield(S(s).gamma,'se') && isfield(S(s).gamma.se,'t0') && isfinite(S(s).gamma.se.t0)
            se_gam(s) = max(0, S(s).gamma.se.t0);
        elseif isfield(S(s).gamma,'ci') && isfield(S(s).gamma.ci,'t0')
            ci = S(s).gamma.ci.t0(:);
            if numel(ci)>=2 && all(isfinite(ci))
                se_gam(s) = (ci(2) - ci(1)) / (2*z);
            else
                se_gam(s) = NaN;
            end
        else
            se_gam(s) = NaN;
        end
    end
end

% tags
Tags = table();
for i=1:numel(tag_fields)
    f = tag_fields{i};
    vals = strings(Ns,1);
    for s=1:Ns
        if isfield(S(s), f) && ~isempty(S(s).(f))
            if isstring(S(s).(f)) || ischar(S(s).(f)) || isscalar(S(s).(f))
                vals(s) = string(S(s).(f));
            else
                try, vals(s) = string(S(s).(f)); catch, vals(s)="MISSING"; end
            end
        else
            vals(s) = "MISSING";
        end
    end
    Tags.(f) = vals;
end

% differences & relative
diff_abs = t0_gam - t0_rate;
SE_diff  = sqrt( max(se_rate,0).^2 + max(se_gam,0).^2 );

rel = nan(Ns,1); rel_SE = nan(Ns,1);
switch lower(string(opts.ratio_mode))
    case {"rate","abs"}  % alias
        denom = max(abs(t0_rate), opts.den_floor);
        rel   = diff_abs ./ denom;
        rel_SE = SE_diff ./ denom;   % ignore denom uncertainty
    case "none"
        % leave rel, rel_SE as NaN; user just uses absolute differences
    otherwise
        % default to 'rate'
        denom = max(abs(t0_rate), opts.den_floor);
        rel   = diff_abs ./ denom;
        rel_SE = SE_diff ./ denom;
end

% CIs
diff_CIlo = diff_abs - z*SE_diff;
diff_CIhi = diff_abs + z*SE_diff;
rel_CIlo  = rel - z*rel_SE;
rel_CIhi  = rel + z*rel_SE;
rel_pct   = 100*rel;
rel_pct_lo= 100*rel_CIlo;
rel_pct_hi= 100*rel_CIhi;

% per-sample table
T_samples = table(Name, t0_rate, se_rate, t0_gam, se_gam, used_fixed_t0, ...
                  diff_abs, SE_diff, diff_CIlo, diff_CIhi, ...
                  rel, rel_SE, rel_CIlo, rel_CIhi, rel_pct, rel_pct_lo, rel_pct_hi);
T_samples = [T_samples, Tags];

% per-condition (exclude reptag for grouping)
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
for k=1:Nc
    rows = (Kidx==k) & isfinite(diff_abs) & isfinite(SE_diff) & (SE_diff > 0);
    baseRow = table();
    for i=1:numel(key_fields)
        f = key_fields{i};
        kk = find(Kidx==k,1,'first');
        baseRow.(f) = Tags.(f)(kk);
    end

    if ~any(rows)
        Row = [ baseRow, table(0, NaN, NaN, NaN, NaN, NaN, NaN, ...
                               'VariableNames', {'Nrep','mean_diff','SE_diff','CIlo_diff','CIhi_diff', ...
                                                 'mean_rel','SE_rel'}) ];
    else
        % absolute difference (inverse-variance)
        w  = 1 ./ (SE_diff(rows).^2);
        w(~isfinite(w)) = 0;
        md = sum(w .* diff_abs(rows)) / sum(w);
        sed = sqrt(1 / sum(w));
        cdlo = md - z*sed; cdhi = md + z*sed;

        % relative (only where rel & rel_SE finite and >0)
        rows_rel = rows & isfinite(rel) & isfinite(rel_SE) & rel_SE > 0;
        if any(rows_rel)
            wr = 1 ./ (rel_SE(rows_rel).^2);
            mr = sum(wr .* rel(rows_rel)) / sum(wr);
            ser = sqrt(1 / sum(wr));
        else
            mr = NaN; ser = NaN;
        end

        Row = [ baseRow, table(sum(rows), md, sed, cdlo, cdhi, mr, ser, ...
                               'VariableNames', {'Nrep','mean_diff','SE_diff','CIlo_diff','CIhi_diff', ...
                                                 'mean_rel','SE_rel'}) ];
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

% write files if requested
if ~isempty(target_folder)
    if ~exist(target_folder,'dir'), mkdir(target_folder); end
    % CSVs
    writetable(T_samples,    fullfile(target_folder, sprintf('%s__samples.csv',    opts.csv_prefix)));
    writetable(T_conditions, fullfile(target_folder, sprintf('%s__conditions.csv', opts.csv_prefix)));
    % Excel (single workbook with two sheets)
    try
        xfile = fullfile(target_folder, opts.excel_name);
        writetable(T_samples,    xfile, 'Sheet','samples');
        writetable(T_conditions, xfile, 'Sheet','conditions');
    catch ME %#ok<NASGU>
        % ignore Excel write issues silently; CSVs are already written
    end
end
end


function v = getfield_soft(S, f, dflt)
if isstruct(S) && isfield(S,f) && ~isempty(S.(f))
    v = S.(f);
else
    v = dflt;
end
end

function G = compute_gof(Ss, modelname)
G = struct('R2',NaN,'AIC',NaN,'BIC',NaN);
if ~isfield(Ss,'bins') || isempty(Ss.bins), return; end
if ~isfield(Ss, modelname) || isempty(Ss.(modelname)) || ~isfield(Ss.(modelname),'fit')
    return;
end
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
ll_null  = sum( y .* log(max(p_bar,eps)) + (n - y) .* log(max(1 - p_bar, eps)) );

if isfinite(ll_model) && isfinite(ll_null) && ll_null ~= 0
    G.R2 = 1 - (ll_model / ll_null);
end
if isfield(Ss.(modelname),'aic'), G.AIC = Ss.(modelname).aic; end
if isfield(Ss.(modelname),'bic'), G.BIC = Ss.(modelname).bic; end
end
