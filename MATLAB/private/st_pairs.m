function pairs = st_pairs(options, names)
% Turn selected struct fields into name/value pairs.
pairs = cell(1, 2*numel(names));
for k = 1:numel(names)
    pairs{2*k-1} = names{k};
    pairs{2*k} = options.(names{k});
end
end
