function out = lookup_calibration(cal_dates, cal_values, varargin)
%LOOKUP_CALIBRATION Match studies to calibration dates, following the R tool.
%   OUT = LOOKUP_CALIBRATION(CAL_DATES, CAL_VALUES, 'filename', FILES)
%   or (..., 'date', DATES). Supply exactly one study input. DATES may also be
%   a SWISSTRACE_CORRECT result. Options: method='last'|'nearest'|'exact',
%   max_gap=[] (days), tz='UTC'. Returns study_date/calibration_factor/cal_date/
%   gap_days in a table. Dates are timezone-free datetime values; missing
%   factors and dates are NaN/NaT. Signed gap = study minus calibration date.
p = inputParser;
p.PartialMatching = false;
p.CaseSensitive = true;
addParameter(p, 'filename', []);
addParameter(p, 'date', []);
addParameter(p, 'method', 'last');
addParameter(p, 'max_gap', []);
addParameter(p, 'tz', 'UTC');
parse(p, varargin{:});
o = p.Results;
has_file = ~ismember('filename', p.UsingDefaults);
has_date = ~ismember('date', p.UsingDefaults);
if has_file == has_date
    error('swisstrace:StudyRequired', 'Supply exactly one of filename or date.');
end
method = st_text(o.method, 'method');
if ~ismember(method, {'last','nearest','exact'})
    error('swisstrace:InvalidMethod', 'method must be last, nearest or exact.');
end
if ~isempty(o.max_gap)
    validateattributes(o.max_gap, {'numeric'}, {'real','scalar','finite','nonnegative'}, mfilename, 'max_gap');
end
if has_file
    files = string(o.filename);
    study_date = NaT(numel(files), 1);
    for k = 1:numel(files), study_date(k) = swisstrace_date(files(k), o.tz); end
    study_date.Format = 'yyyy-MM-dd';
else
    study_date = st_dates(o.date);
end
cal_dates = st_dates(cal_dates);
if ~isnumeric(cal_values), cal_values = str2double(string(cal_values)); end
cal_values = cal_values(:);
if numel(cal_dates) ~= numel(cal_values)
    error('swisstrace:CalibrationLength', 'cal_dates and cal_values must have the same length.');
end
if any(isnat(cal_dates))
    error('swisstrace:InvalidDate', 'cal_dates contains unparseable dates.');
end
[cal_dates, order] = sort(cal_dates);
cal_values = cal_values(order);
n = numel(study_date);
calibration_factor = NaN(n, 1);
cal_date = NaT(n, 1);
cal_date.Format = 'yyyy-MM-dd';
gap_days = NaN(n, 1);
for k = 1:n
    if isnat(study_date(k)) || isempty(cal_dates), continue; end
    switch method
        case 'last'
            hit = find(cal_dates <= study_date(k), 1, 'last');
        case 'exact'
            hit = find(cal_dates == study_date(k), 1, 'last');
        case 'nearest'
            [~, hit] = min(abs(days(study_date(k) - cal_dates)));
    end
    if isempty(hit), continue; end
    cal_date(k) = cal_dates(hit);
    gap_days(k) = days(study_date(k) - cal_date(k));
    if isempty(o.max_gap) || abs(gap_days(k)) <= o.max_gap
        calibration_factor(k) = cal_values(hit);
    end
end
if any(isnan(calibration_factor))
    warning('swisstrace:MissingCalibration', '%d of %d studies had no calibration factor within range.', sum(isnan(calibration_factor)), n);
end
out = table(study_date, calibration_factor, cal_date, gap_days);
end
