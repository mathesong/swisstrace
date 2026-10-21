function st_write_curve(path, header, tac, missing_value)
fid = fopen(path, 'wt');
if fid < 0, error('swisstrace:WriteFailed', 'Cannot write %s.', path); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', header);
for k = 1:height(tac)
    if isfinite(tac.activity(k))
        fprintf(fid, '%.17g\t%.17g\n', tac.time(k), tac.activity(k));
    else
        fprintf(fid, '%.17g\t%s\n', tac.time(k), missing_value);
    end
end
end
