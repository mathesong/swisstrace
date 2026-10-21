function o = st_options(kind, varargin)
% Shared defaults; reject misspelled or partial option names.
o = struct('pet_start', [], 'half_life', [], 'lead', 20, ...
    'frame_scheme', [], 'zero_first_frame', true, 'baseline_k', 3, ...
    'min_run', 5, 'baseline_init', 10, 'baseline_min', 20, 'tz', 'UTC');
if any(strcmp(kind, {'process', 'batch'}))
    o.output_folder = [];
    o.sub = [];
    o.ses = [];
    o.bids_dir = [];
    o.recording = 'autosampler';
    o.overwrite = true;
end
if strcmp(kind, 'batch')
    o.cal_dates = [];
    o.cal_values = [];
    o.cal_method = 'last';
    o.cal_max_gap = [];
end
p = inputParser;
p.PartialMatching = false;
p.CaseSensitive = true;
names = fieldnames(o);
for k = 1:numel(names)
    addParameter(p, names{k}, o.(names{k}));
end
parse(p, varargin{:});
o = p.Results;
end
