function status_table = swisstrace_convert_batch(manifest, varargin)
%SWISSTRACE_CONVERT_BATCH Convert CSV/TSV/Excel manifests or a MATLAB table.
%   STATUS = SWISSTRACE_CONVERT_BATCH(MANIFEST, Name, Value, ...) accepts
%   output_folder, cal_dates, cal_values, bids_dir, cal_method='last',
%   cal_max_gap=[], and all SWISSTRACE_PROCESS options. Manifest columns:
%   filename, isotope (or half_life), pet_start, calibration_factor, sub, ses.
%   Column names are case-insensitive. Relative files resolve from the
%   manifest folder (pwd for table input). File labels retain leading zeros.
%   Missing factors use lookup_calibration and are logged. Each bad row is
%   reported without aborting other rows, including failed date lookups.
o = st_options('batch', varargin{:});
if ~ismember(string(o.cal_method), ["last","nearest","exact"])
    error('swisstrace:InvalidMethod', 'cal_method must be last, nearest or exact.');
end
if istable(manifest)
    man = manifest;
    manifest_dir = pwd;
else
    path = st_text(manifest, 'manifest');
    if ~isfile(path), error('swisstrace:FileNotFound', 'Manifest not found: %s', path); end
    [~, attr] = fileattrib(path);
    manifest_dir = fileparts(attr.Name);
    [~, ~, ext] = fileparts(path);
    switch lower(ext)
        case {'.xlsx','.xls'}
            opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');
        case {'.csv','.tsv','.tab','.txt'}
            delimiter = '\t';
            if strcmpi(ext, '.csv'), delimiter = ','; end
            opts = detectImportOptions(path, 'FileType', 'text', ...
                'Delimiter', delimiter, 'VariableNamingRule', 'preserve');
        otherwise
            error('swisstrace:ManifestType', 'Unsupported manifest type; use CSV, TSV, XLSX or XLS.');
    end
    opts = setvartype(opts, opts.VariableNames, 'string');
    man = readtable(path, opts);
end
names = lower(strtrim(string(man.Properties.VariableNames)));
if numel(unique(names)) ~= numel(names)
    error('swisstrace:ManifestColumns', 'Duplicate manifest column names after normalization.');
end
man.Properties.VariableNames = cellstr(names);
if ~ismember('filename', names)
    error('swisstrace:ManifestColumns', 'Manifest is missing required column: filename.');
end
root = o.output_folder;
if isempty(root), root = manifest_dir; end
root = st_text(root, 'output_folder');
n = height(man);
filename = strings(n, 1);
isotope = strings(n, 1);
calibration_factor = NaN(n, 1);
for k = 1:n
    filename(k) = text_cell(man, 'filename', k);
    isotope(k) = text_cell(man, 'isotope', k);
    calibration_factor(k) = number_cell(man, 'calibration_factor', k);
end
need = isnan(calibration_factor);
if any(need) && (isempty(o.cal_dates) || isempty(o.cal_values))
    error('swisstrace:CalibrationRequired', 'Missing calibration_factor; supply cal_dates and cal_values.');
end
if ~isfolder(root), mkdir(root); end
status = repmat("ok", n, 1);
message = strings(n, 1);
log_date = NaT(n, 1);
log_date.Format = 'yyyy-MM-dd';
log_gap = NaN(n, 1);
for k = 1:n
    try
        if ismissing(filename(k)) || strlength(filename(k)) == 0
            error('swisstrace:FileNotFound', 'filename is missing.');
        end
        path = char(filename(k));
        % Absolute Unix, drive-letter and UNC paths; other paths are relative.
        if isempty(regexp(path, '^(/|[A-Za-z]:[\\/]|\\\\)', 'once'))
            path = fullfile(manifest_dir, path);
        end
        if need(k)
            match = lookup_calibration(o.cal_dates, o.cal_values, 'filename', path, ...
                'method', o.cal_method, 'max_gap', o.cal_max_gap, 'tz', o.tz);
            calibration_factor(k) = match.calibration_factor;
            log_date(k) = match.study_date;
            log_gap(k) = match.gap_days;
        end
        if ~isfinite(calibration_factor(k))
            error('swisstrace:MissingCalibration', 'No calibration factor found within range.');
        end
        row_options = o;
        row_options.output_folder = root;
        for name = {'sub','ses','pet_start'}
            if ismember(name{1}, names)
                row_options.(name{1}) = cell_value(man, name{1}, k);
            end
        end
        value = row_options.pet_start;
        if ischar(value) || isstring(value)
            number = str2double(value);
            if isfinite(number), row_options.pet_start = number; end
        end
        half_life = number_cell(man, 'half_life', k);
        if ~isnan(half_life), row_options.half_life = half_life; end
        if isempty(row_options.sub), row_options.bids_dir = []; end
        pairs = st_pairs(row_options, fieldnames(st_options('process')));
        iso = isotope(k);
        if ismissing(iso) || strlength(iso) == 0, iso = []; end
        swisstrace_process(path, calibration_factor(k), iso, pairs{:});
    catch err
        status(k) = "error";
        message(k) = string(err.message);
    end
end
if any(need)
    date = string(log_date(need), 'yyyy-MM-dd');
    date(ismissing(date)) = "NA";
    log = table(filename(need), date, calibration_factor(need), log_gap(need), ...
        'VariableNames', {'filename','date','calibration_factor','gap_days'});
    writetable(log, fullfile(root, 'calibration_factors.tsv'), 'FileType', 'text', 'Delimiter', '\t');
end
status_table = table(filename, isotope, calibration_factor, status, message);
end

function value = cell_value(man, name, row)
value = [];
if ~ismember(name, man.Properties.VariableNames), return; end
column = man.(name);
value = column(row, :);
if iscell(value), value = value{1}; end
if isempty(value) || all(ismissing(value)), value = []; return; end
if ischar(value) || isstring(value)
    value = strtrim(string(value));
    if strlength(value) == 0 || value == "NA", value = []; end
end
end

function value = text_cell(man, name, row)
v = cell_value(man, name, row);
value = "";
if ~isempty(v), value = string(v); end
end

function value = number_cell(man, name, row)
v = cell_value(man, name, row);
value = NaN;
if isempty(v), return; end
if isnumeric(v) && isscalar(v), value = double(v); else, value = str2double(string(v)); end
end
