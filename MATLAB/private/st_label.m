function label = st_label(value, entity)
label = regexprep(st_text(string(value), entity), ['^' entity '-'], '');
label = regexprep(label, '[^A-Za-z0-9]', '');
if isempty(label)
    error('swisstrace:InvalidLabel', '%s must contain an alphanumeric BIDS label.', entity);
end
end
