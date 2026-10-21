function lines = st_read_lines(file, limit)
file = st_text(file, 'file');
fid = fopen(file, 'rt');
if fid < 0
    error('swisstrace:FileNotFound', 'Cannot open file: %s', file);
end
cleanup = onCleanup(@() fclose(fid));
lines = cell(0, 1);
while numel(lines) < limit
    line = fgetl(fid);
    if ~ischar(line), break; end
    lines{end+1, 1} = strtrim(line); %#ok<AGROW>
end
end
