function value = swisstrace_date(file, tz)
%SWISSTRACE_DATE Read acquisition date from the first timestamp in 50 lines.
%   DATE = SWISSTRACE_DATE(FILE, TZ) returns a timezone-free datetime at midnight.
%   TZ defaults to 'UTC'. Like R, zoned timestamps are converted to a UTC date.
if nargin < 2, tz = 'UTC'; end
lines = st_read_lines(file, 50);
for k = 1:numel(lines)
    parts = regexp(lines{k}, '\s+', 'split');
    if numel(parts) < 6, continue; end
    v = str2double(parts(1:6));
    if all(isfinite(v))
        value = st_dates(st_timestamp(v, tz));
        return;
    end
end
error('swisstrace:InvalidTimestamp', 'Could not read a timestamp from %s.', file);
end
