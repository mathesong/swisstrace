function res = swisstrace_correct(file, calibration_factor, isotope, varargin)
%SWISSTRACE_CORRECT Correct/calibrate a raw Swisstrace twilite recording.
%   RES = SWISSTRACE_CORRECT(FILE, CALIBRATION_FACTOR, ISOTOPE, Name, Value, ...)
%   returns a struct containing tac/raw tables and correction metadata.
%   Use ISOTOPE=[] with 'half_life' to supply a half-life directly in seconds.
%
%   Options (R defaults): pet_start=[], half_life=[], lead=20, frame_scheme=[],
%   zero_first_frame=true, baseline_k=3, min_run=5, baseline_init=10,
%   baseline_min=20, tz='UTC'. pet_start accepts elapsed seconds, datetime, or
%   'HH:MM[:SS[.sss]]' on the acquisition date. [] detects the bolus onset.
%   frame_scheme is a table with width/end columns, or an N-by-2 [width end]
%   matrix. No framing is performed by default. Empty frames contain NaN.
%
%   res = swisstrace_correct('raw.crv', .425, 'F18', 'pet_start', '11:16:30');
if nargin < 3, isotope = []; end
if nargin < 2
    error('swisstrace:CalibrationRequired', 'calibration_factor is required.');
end
o = st_options('correct', varargin{:});
assert_raw_crv(file);
validateattributes(calibration_factor, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'calibration_factor');
for name = {'lead','baseline_k','baseline_min'}
    validateattributes(o.(name{1}), {'numeric'}, {'real','scalar','finite','nonnegative'}, mfilename, name{1});
end
validateattributes(o.baseline_init, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'baseline_init');
validateattributes(o.min_run, {'numeric'}, {'real','scalar','finite','positive','integer'}, mfilename, 'min_run');
validateattributes(o.zero_first_frame, {'logical','numeric'}, {'scalar','binary'}, mfilename, 'zero_first_frame');
names = {'F18','C11','N13','O15','Ga68','Cu62','Zr89','I124','Rb82'};
lives = [6586.2,1223.4,597.9,122.24,4057.74,584.4,282240,360806.4,76.4];
if isempty(o.half_life)
    if isempty(isotope)
        error('swisstrace:IsotopeRequired', 'isotope is required unless half_life is supplied.');
    end
    isotope = st_text(isotope, 'isotope');
    key = regexprep(isotope, '[^A-Za-z0-9]', '');
    index = find(strcmpi(key, names), 1);
    if isempty(index)
        error('swisstrace:UnknownIsotope', 'Unknown isotope %s; supply half_life.', isotope);
    end
    half_life = lives(index);
else
    half_life = o.half_life;
end
validateattributes(half_life, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'half_life');
lambda = log(2) / half_life;
lines = st_read_lines(file, Inf);
m = zeros(numel(lines), 9);
n = 0;
for k = 1:numel(lines)
    parts = regexp(lines{k}, '\s+', 'split');
    if numel(parts) < 9, continue; end
    n = n + 1;
    m(n, :) = str2double(parts(1:9));
end
m = m(1:n, :);
if n == 0 || any(~isfinite(m(:)))
    error('swisstrace:InvalidRawData', 'No valid data rows, or nonfinite timestamps/counts.');
end
ts = st_timestamp(m(:, 1:6), o.tz);
secs = seconds(ts - ts(1));
if any(diff(secs) < 0)
    error('swisstrace:UnsortedTime', 'Raw timestamps must be chronological.');
end
coinc = m(:, 7);
detected = isempty(o.pet_start);
lead_used = NaN;
if detected
    bg = find(secs < o.baseline_init);
    if numel(bg) < 3, bg = (1:min(3, n))'; end
    onset = [];
    for k = (max(bg) + 1):n
        baseline = coinc(bg);
        % R stats::mad uses 1.4826 times the median absolute deviation.
        spread = 1.4826 * median(abs(baseline - median(baseline)));
        if ~isfinite(spread) || spread == 0, spread = std(baseline, 0); end
        if ~isfinite(spread) || spread == 0, spread = 1; end
        threshold = median(baseline) + o.baseline_k * spread;
        last = k + o.min_run - 1;
        if last <= n && all(coinc(k:last) > threshold)
            onset = k;
            break;
        end
        if coinc(k) <= threshold, bg(end+1, 1) = k; end %#ok<AGROW>
    end
    if isempty(onset)
        warning('swisstrace:NoRise', 'No sustained rise above background; using recording start as time 0.');
        t0 = 0;
    else
        onset_sec = secs(onset);
        if onset_sec < o.baseline_min
            warning('swisstrace:ShortBaseline', 'Rise after %.1f s; background may be unreliable.', onset_sec);
        end
        lead_used = min(o.lead, onset_sec);
        if lead_used < o.lead
            warning('swisstrace:ReducedLead', 'Only %.1f s before rise; using a %.1f s lead.', onset_sec, lead_used);
        end
        t0 = onset_sec - lead_used;
    end
else
    start = o.pet_start;
    if ischar(start) || isstring(start)
        clock = strtrim(st_text(start, 'pet_start'));
        if isempty(regexp(clock, '^\d{1,2}:\d{2}(:\d{2}(\.\d+)?)?$', 'once'))
            error('swisstrace:InvalidClock', 'pet_start must be a valid HH:MM[:SS] clock time.');
        end
        v = str2double(strsplit(clock, ':'));
        if numel(v) == 2, v(3) = 0; end
        if v(1) > 23 || v(2) > 59 || v(3) >= 60
            error('swisstrace:InvalidClock', 'pet_start must be a valid clock time.');
        end
        start = datetime(year(ts(1)), month(ts(1)), day(ts(1)), v(1), v(2), v(3), 'TimeZone', o.tz);
    end
    if isdatetime(start)
        if ~isscalar(start) || isnat(start)
            error('swisstrace:InvalidStart', 'pet_start must be a valid scalar datetime.');
        end
        if isempty(start.TimeZone), start.TimeZone = o.tz; end
        t0 = seconds(start - ts(1));
    else
        validateattributes(start, {'numeric'}, {'real','scalar','finite'}, mfilename, 'pet_start');
        t0 = start;
    end
    if t0 < 0
        warning('swisstrace:ClampedStart', 'Time 0 precedes recording start; clamping to 0.');
        t0 = 0;
    end
    bg = find(secs < t0);
end
if t0 > secs(end)
    error('swisstrace:StartAfterEnd', 'pet_start is after the end of the recording.');
end
if isempty(bg)
    background = 0;
    warning('swisstrace:NoBackground', 'No samples before supplied time 0; background set to 0.');
else
    background = mean(coinc(bg));
end
tau = secs - t0;
keep = tau >= 0;
time = tau(keep);
activity = (coinc(keep) - background) .* exp(lambda * time) * calibration_factor;
if isempty(o.frame_scheme)
    if numel(time) >= 2
        mids = (time(1:end-1) + time(2:end)) / 2;
        frame_start = [0; mids];
        frame_end = [mids; 2*time(end) - mids(end)];
    else
        frame_start = zeros(size(time));
        frame_end = time;
    end
else
    scheme = o.frame_scheme;
    if istable(scheme)
        if ~all(ismember({'width','end'}, scheme.Properties.VariableNames))
            error('swisstrace:InvalidScheme', 'frame_scheme needs width and end columns.');
        end
        scheme = [scheme.width, scheme.('end')];
    end
    validateattributes(scheme, {'numeric'}, {'real','2d','ncols',2,'finite','positive'}, mfilename, 'frame_scheme');
    edges = 0;
    previous = 0;
    for k = 1:size(scheme, 1)
        width = scheme(k, 1); stop = scheme(k, 2);
        if stop - previous < width
            error('swisstrace:InvalidScheme', 'Frame ends must increase by at least their width.');
        end
        edges = [edges; (previous + width:width:stop)']; %#ok<AGROW>
        previous = stop;
    end
    edges = unique(edges);
    nframes = min(sum(edges <= time(end)), numel(edges) - 1);
    values = NaN(nframes, 1);
    for k = 1:nframes
        selected = time >= edges(k) & time < edges(k+1);
        if any(selected), values(k) = mean(activity(selected)); end
    end
    frame_start = edges(1:nframes);
    frame_end = edges(2:nframes+1);
    time = (frame_start + frame_end) / 2;
    activity = values;
end
if o.zero_first_frame && ~isempty(activity), activity(1) = 0; end
frame_dur = frame_end - frame_start;
res = struct;
res.tac = table(time, activity, frame_start, frame_end, frame_dur);
res.calibration_factor = calibration_factor;
res.background = background;
res.isotope = isotope;
res.half_life = half_life;
res.lambda = lambda;
res.date = st_dates(ts(1));
res.start_time = char(string(ts(1), 'HH:mm:ss.SSS'));
res.acq_start = ts(1);
res.acq_end = ts(end);
res.acq_duration = secs(end);
res.pet_start = ts(1) + seconds(t0);
res.t0_seconds = t0;
res.t0_detected = detected;
res.lead = lead_used;
res.n_raw = n;
res.n_background = numel(bg);
[~, attr] = fileattrib(file);
res.file = attr.Name;
res.frame_scheme = o.frame_scheme;
res.raw = table(secs, coinc, m(:, 8), m(:, 9), ...
    'VariableNames', {'time','coincidence','singles1','singles2'});
end
