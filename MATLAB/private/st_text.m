function value = st_text(value, name)
if ~(ischar(value) && isrow(value)) && ...
        ~(isstring(value) && isscalar(value) && ~ismissing(value))
    error('swisstrace:InvalidText', '%s must be a text scalar.', name);
end
value = char(value);
end
