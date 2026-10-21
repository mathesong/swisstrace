function valid = assert_raw_crv(file)
%ASSERT_RAW_CRV Reject empty or corrected/processed twilite recordings.
%   ASSERT_RAW_CRV(FILE) returns true for a numeric raw header, otherwise errors.
%   Full row validation is performed by SWISSTRACE_CORRECT.
lines = st_read_lines(file, 5);
lines = lines(~cellfun('isempty', lines));
if isempty(lines)
    error('swisstrace:EmptyFile', 'File is empty: %s', file);
end
first = lines{1};
parts = regexp(first, '\s+', 'split');
if ~isfinite(str2double(parts{1})) || ...
        ~isempty(regexpi(first, 'corrected|kBq/cc|value\[|time\[seconds\]', 'once'))
    error('swisstrace:NotRaw', ...
        ['File looks like a corrected/processed .crv, not a raw twilite recording. ' ...
         'Expected YYYY M D H M S coincidences singles1 singles2.']);
end
valid = true;
end
