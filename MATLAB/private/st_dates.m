function dates = st_dates(value)
% Calendar dates are timezone-free datetimes. Match R's UTC POSIXct -> Date.
if isstruct(value) && isfield(value, 'date'), value = value.date; end
if isdatetime(value)
    if ~isempty(value.TimeZone), value.TimeZone = 'UTC'; end
    dates = dateshift(value(:), 'start', 'day');
    dates.TimeZone = '';
elseif isempty(value)
    dates = NaT(0, 1);
elseif ischar(value) || isstring(value) || iscellstr(value)
    value = string(value);
    dates = NaT(numel(value), 1);
    for k = 1:numel(value)
        if ismissing(value(k)), continue; end
        tokens = regexp(char(strtrim(value(k))), '^(\d{4})[-/](\d{1,2})[-/](\d{1,2})$', 'tokens', 'once');
        if isempty(tokens)
            error('swisstrace:InvalidDate', 'Use YYYY-MM-DD or YYYY/MM/DD for study dates.');
        end
        v = str2double(tokens);
        ts = st_timestamp([v, 0, 0, 0], 'UTC');
        ts.TimeZone = '';
        dates(k) = ts;
    end
else
    error('swisstrace:InvalidDate', 'Dates must be datetime values or date strings.');
end
dates.Format = 'yyyy-MM-dd';
end
