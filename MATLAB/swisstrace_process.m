function written = swisstrace_process(filename, calibration_factor, isotope, varargin)
%SWISSTRACE_PROCESS Correct a raw CRV and write PMOD, PNG and optional BIDS.
%   PATHS = SWISSTRACE_PROCESS(FILE, FACTOR, ISOTOPE, Name, Value, ...)
%   accepts all SWISSTRACE_CORRECT options plus output_folder=[], sub=[],
%   ses=[], bids_dir=[], recording='autosampler', overwrite=true.
%   Subject/session labels are strings (e.g. '01'). With bids_dir, a matching
%   PET image supplies the output directory/entities; multiple matches warn
%   and use the first sorted path, as in R. Returns a struct of output paths.
if nargin < 3, isotope = []; end
o = st_options('process', varargin{:});
assert_raw_crv(filename);
if isempty(o.sub) && (~isempty(o.ses) || ~isempty(o.bids_dir))
    error('swisstrace:SubjectRequired', 'ses or bids_dir was provided without sub.');
end
validateattributes(o.overwrite, {'logical','numeric'}, {'scalar','binary'}, mfilename, 'overwrite');
if ~isempty(o.sub)
    subject = st_label(o.sub, 'sub');
    session = '';
    if ~isempty(o.ses), session = st_label(o.ses, 'ses'); end
    recording = st_label(o.recording, 'recording');
end
names = fieldnames(st_options('correct'));
pairs = st_pairs(o, names);
res = swisstrace_correct(filename, calibration_factor, isotope, pairs{:});
[source_folder, base, ~] = fileparts(filename);
root = o.output_folder;
if isempty(root), root = source_folder; end
if isempty(root), root = pwd; end
root = st_text(root, 'output_folder');
written = struct('corrected_crv', fullfile(root, 'Corrected_PMOD', [base '_corrected.crv']), ...
    'plot_png', fullfile(root, 'Plots', [base '_corrected.png']));
if ~isempty(o.sub)
    bids = o.bids_dir;
    if isempty(bids), bids = fullfile(root, 'BIDS'); end
    pet_dir = fullfile(bids, ['sub-' subject]);
    entities = ['sub-' subject];
    if ~isempty(session)
        pet_dir = fullfile(pet_dir, ['ses-' session]);
        entities = [entities '_ses-' session];
    end
    pet_dir = fullfile(pet_dir, 'pet');
    if ~isempty(o.bids_dir)
        files = dir(fullfile(bids, '**', '*_pet.nii*'));
        hits = strings(0, 1);
        pattern = ['^' entities '(_.*)?_pet\.nii(\.gz)?$'];
        for k = 1:numel(files)
            if ~files(k).isdir && ~isempty(regexp(files(k).name, pattern, 'once'))
                hits(end+1, 1) = string(fullfile(files(k).folder, files(k).name)); %#ok<AGROW>
            end
        end
        hits = sort(hits);
        if ~isempty(hits)
            if numel(hits) > 1
                warning('swisstrace:MultiplePET', '%d PET images matched; using first: %s', numel(hits), hits(1));
            end
            [pet_dir, name, ext] = fileparts(char(hits(1)));
            entities = regexprep([name ext], '_pet\.nii(\.gz)?$', '');
        end
    end
    stem = [entities '_recording-' recording '_blood'];
    written.bids_tsv = fullfile(pet_dir, [stem '.tsv']);
    written.bids_json = fullfile(pet_dir, [stem '.json']);
end
names = fieldnames(written);
for k = 1:numel(names)
    kind = names{k}; path = written.(kind);
    folder = fileparts(path);
    if ~isfolder(folder), mkdir(folder); end
    if isfile(path) && ~o.overwrite, continue; end
    switch kind
        case 'corrected_crv'
            header = sprintf('Corrected_&_calibrated_[kBq/cc]_>___time[seconds]\tvalue[kBq/cc]');
            st_write_curve(path, header, res.tac, 'NA');
        case 'plot_png'
            fig = swisstrace_qc(res, 'visible', 'off');
            cleanup = onCleanup(@() close(fig));
            print(fig, path, '-dpng', '-r110');
            clear cleanup;
        case 'bids_tsv'
            st_write_curve(path, sprintf('time\twhole_blood_radioactivity'), res.tac, 'n/a');
        case 'bids_json'
            meta = struct('PlasmaAvail', false, 'WholeBloodAvail', true, ...
                'MetaboliteAvail', false, 'DispersionCorrected', false);
            meta.time = struct('Description', ...
                'Time relative to time zero (PET scan start); native sample time or frame midpoint.', 'Units', 's');
            meta.whole_blood_radioactivity = struct('Description', ...
                'Whole blood radioactivity concentration: background-subtracted, decay-corrected to time zero, and calibrated.', ...
                'Units', 'kBq/mL');
            fid = fopen(path, 'wt');
            if fid < 0, error('swisstrace:WriteFailed', 'Cannot write %s.', path); end
            cleanup = onCleanup(@() fclose(fid));
            fprintf(fid, '%s\n', jsonencode(meta, 'PrettyPrint', true));
            clear cleanup;
    end
end
end
