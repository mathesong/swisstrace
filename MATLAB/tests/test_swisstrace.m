function tests = test_swisstrace
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.root = root;
testCase.TestData.old_path = path;
addpath(root);
end

function teardownOnce(testCase)
path(testCase.TestData.old_path);
end

function setup(testCase)
testCase.TestData.folder = tempname;
mkdir(testCase.TestData.folder);
end

function teardown(testCase)
rmdir(testCase.TestData.folder, 's');
end

function testExactArithmeticAndMetadata(testCase)
f = exact_file(testCase);
r = swisstrace_correct(f, .5, 'F18', 'pet_start', 40);
verifyEqual(testCase, r.background, 30);
verifyEqual(testCase, r.n_background, 40);
verifyEqual(testCase, r.n_raw, 301);
verifyEqual(testCase, r.t0_seconds, 40);
verifyFalse(testCase, r.t0_detected);
verifyEqual(testCase, r.tac.activity(1), 0);
verifyEqual(testCase, r.tac.activity(2:end), repmat(50, 260, 1), 'AbsTol', 1e-10);
verifyEqual(testCase, r.tac.Properties.VariableNames, {'time','activity','frame_start','frame_end','frame_dur'});
verifyEqual(testCase, r.raw.Properties.VariableNames, {'time','coincidence','singles1','singles2'});
verifyEqual(testCase, r.tac.frame_start(2:end), r.tac.frame_end(1:end-1));
verifyEqual(testCase, r.date, datetime(2026,1,1));
verifyEqual(testCase, swisstrace_date(f), r.date);
end

function testStartFormatsAndHalfLifeOverride(testCase)
f = exact_file(testCase);
starts = {40, '00:00:40', datetime(2026,1,1,0,0,40,'TimeZone','UTC'), datetime(2026,1,1,0,0,40)};
for k = 1:numel(starts)
    r = swisstrace_correct(f, .5, [], 'half_life', 6586.2, ...
        'pet_start', starts{k}, 'zero_first_frame', false);
    verifyEqual(testCase, r.tac.activity, repmat(50,261,1), 'AbsTol', 1e-10);
    verifyEqual(testCase, r.t0_seconds, 40);
end
end

function testIsotopeTable(testCase)
f = exact_file(testCase);
names = {'f-18','c11','n13','o15','ga68','cu62','zr89','i124','rb82'};
expected = [6586.2,1223.4,597.9,122.24,4057.74,584.4,282240,360806.4,76.4];
for k = 1:numel(names)
    r = swisstrace_correct(f, .5, names{k}, 'pet_start', 40);
    verifyEqual(testCase, r.half_life, expected(k));
end
end

function testAutoOnsetAndSpikeExclusion(testCase)
f = exact_file(testCase);
r = swisstrace_correct(f, .5, 'F18');
verifyEqual(testCase, r.t0_seconds, 20);
verifyEqual(testCase, r.lead, 20);
verifyEqual(testCase, r.background, 30);
verifyTrue(testCase, r.t0_detected);
counts = [repmat(30,40,1); repmat(130,30,1)];
counts(21) = 1000;
f = write_raw(testCase.TestData.folder, (0:69)', counts);
r = swisstrace_correct(f, .5, 'F18');
verifyEqual(testCase, r.background, 30);
verifyEqual(testCase, r.n_background, 39);
verifyEqual(testCase, r.t0_seconds, 20);
end

function testShortFlatRecording(testCase)
f = write_raw(testCase.TestData.folder, [0;1], [30;30]);
verifyWarning(testCase, @() swisstrace_correct(f,.5,'F18'), 'swisstrace:NoRise');
end

function testNativeIrregularAndFractionalTimes(testCase)
t = [0;1.25;2.5;4.75];
f = write_raw(testCase.TestData.folder, t, [30;130;120;110]);
r = swisstrace_correct(f, .5, 'F18', 'pet_start', 1, 'zero_first_frame', false);
verifyEqual(testCase, r.tac.time, [.25;1.5;3.75], 'AbsTol', 1e-9);
verifyEqual(testCase, r.tac.frame_start, [0;.875;2.625], 'AbsTol', 1e-9);
verifyEqual(testCase, r.tac.frame_end, [.875;2.625;4.875], 'AbsTol', 1e-9);
end

function testFramesAndEmptyBins(testCase)
t = [0;1;2;4;6];
counts = 30 + [0;1;3;5;7] .* exp(-log(2)/6586.2*(t-1));
f = write_raw(testCase.TestData.folder,t,counts);
scheme = table([1;2],[3;7],'VariableNames',{'width','end'});
r = swisstrace_correct(f,1,'F18','pet_start',1,'zero_first_frame',false,'frame_scheme',scheme);
verifyEqual(testCase, r.tac.time, [.5;1.5;2.5;4;6]);
verifyTrue(testCase, isnan(r.tac.activity(3)));
verifyEqual(testCase, r.tac.activity([1,2,4,5]), [1;3;5;7], 'AbsTol', 1e-10);
% The final right endpoint must be excluded from the last frame, as in R.
r = swisstrace_correct(f,1,'F18','pet_start',1,'zero_first_frame',false,'frame_scheme',[1 5]);
verifyEqual(testCase, height(r.tac), 5);
verifyTrue(testCase, isnan(r.tac.activity(end)));
end

function testValidation(testCase)
f = exact_file(testCase);
verifyError(testCase, @() swisstrace_correct(f,.5), 'swisstrace:IsotopeRequired');
verifyError(testCase, @() swisstrace_correct(f,.5,'XX'), 'swisstrace:UnknownIsotope');
verifyError(testCase, @() swisstrace_correct(f,.5,'F18','pet_start','25:00'), 'swisstrace:InvalidClock');
verifyError(testCase, @() swisstrace_correct(f,.5,'F18','pet_start',9999), 'swisstrace:StartAfterEnd');
verifyError(testCase, @() swisstrace_correct(f,.5,'F18','frame_scheme',[5 3]), 'swisstrace:InvalidScheme');
verifyError(testCase, @() swisstrace_process(f,.5,'F18','ses','01'), 'swisstrace:SubjectRequired');
fid = fopen(f, 'wt'); fprintf(fid, 'Corrected_&_calibrated_[kBq/cc]\n0.5 0\n'); fclose(fid);
verifyError(testCase, @() assert_raw_crv(f), 'swisstrace:NotRaw');
end

function testCalibrationMethodsAndTies(testCase)
dates = ["2026-01-03","2026-01-01","2026-01-01"];
values = [.6,.4,.5];
a = lookup_calibration(dates,values,'date','2026/01/02');
b = lookup_calibration(dates,values,'date','2026-01-02','method','nearest');
c = lookup_calibration(dates,values,'filename',exact_file(testCase),'method','exact');
verifyEqual(testCase, a.calibration_factor, .5);
verifyEqual(testCase, b.calibration_factor, .4);
verifyEqual(testCase, c.calibration_factor, .5);
verifyEqual(testCase, a.gap_days, 1);
verifyError(testCase, @() lookup_calibration(dates,values), 'swisstrace:StudyRequired');
r = swisstrace_correct(exact_file(testCase),.5,'F18');
a = lookup_calibration(dates,values,'date',r);
verifyEqual(testCase, a.calibration_factor, .5);
end

function testCalibrationMissingAndMaxGap(testCase)
old = warning('off','swisstrace:MissingCalibration');
cleanup = onCleanup(@() warning(old));
a = lookup_calibration("2026-01-03",.6,'date',["2025-01-01","2026-02-01"],'max_gap',2);
verifyTrue(testCase, all(isnan(a.calibration_factor)));
verifyEqual(testCase, a.gap_days(2), 29);
verifyEqual(testCase, a.cal_date(2), datetime(2026,1,3));
b = lookup_calibration([],[],'date','2026-01-01');
verifyTrue(testCase, isnan(b.calibration_factor));
end

function testDatesAndMidnight(testCase)
f = write_raw(testCase.TestData.folder,[0;1;2;3],[30;30;130;130], datetime(2026,1,1,23,59,59,'TimeZone','UTC'));
r = swisstrace_correct(f,.5,'F18','pet_start',datetime(2026,1,2,0,0,1,'TimeZone','UTC'));
verifyEqual(testCase, r.t0_seconds, 2);
verifyEqual(testCase, r.acq_duration, 3);
% The R date conversion uses UTC even when recording wall times use another zone.
verifyEqual(testCase, swisstrace_date(exact_file(testCase),'Europe/Copenhagen'),datetime(2025,12,31));
end

function testExportsAndOverwrite(testCase)
f = exact_file(testCase);
p = swisstrace_process(f,.5,'F18','pet_start',40,'sub','sub-01','ses','ses-02');
verifyTrue(testCase, all(structfun(@isfile,p)));
a = readmatrix(p.corrected_crv,'FileType','text','NumHeaderLines',1);
verifyEqual(testCase, a(2:end,2), repmat(50,260,1), 'AbsTol', 1e-10);
meta = jsondecode(fileread(p.bids_json));
verifyTrue(testCase, meta.WholeBloodAvail);
verifyFalse(testCase, meta.PlasmaAvail);
verifyFalse(testCase, meta.MetaboliteAvail);
verifyFalse(testCase, meta.DispersionCorrected);
verifyEqual(testCase, meta.whole_blood_radioactivity.Units, 'kBq/mL');
image = imread(p.plot_png);
verifyGreaterThan(testCase, numel(image), 1000);
before = fileread(p.corrected_crv);
swisstrace_process(f,.7,'F18','pet_start',40,'sub','01','ses','02','overwrite',false);
verifyEqual(testCase, fileread(p.corrected_crv), before);
fig = swisstrace_qc(swisstrace_correct(f,.5,'F18'),'visible','off');
cleanup = onCleanup(@() close(fig));
verifyEqual(testCase, numel(findall(fig,'Type','axes')), 2);
end

function testBidsImageEntitiesAndSubjectBoundary(testCase)
f = exact_file(testCase);
bids = fullfile(testCase.TestData.folder,'dataset');
folder = fullfile(bids,'sub-10','pet'); mkdir(folder);
touch(fullfile(folder,'sub-10_pet.nii.gz'));
p = swisstrace_process(f,.5,'F18','sub','1','bids_dir',bids);
verifyEqual(testCase, fileparts(p.bids_tsv), fullfile(bids,'sub-1','pet'));
folder = fullfile(bids,'sub-01','ses-02','pet'); mkdir(folder);
touch(fullfile(folder,'sub-01_ses-02_trc-FDG_run-1_pet.nii.gz'));
p = swisstrace_process(f,.5,'F18','sub','01','ses','02','bids_dir',bids);
verifyEqual(testCase, char(p.bids_tsv), fullfile(folder,'sub-01_ses-02_trc-FDG_run-1_recording-autosampler_blood.tsv'));
end

function testBatchFormatsAndFailures(testCase)
f = exact_file(testCase);
[~,base,ext] = fileparts(f);
filename = [string([base ext]);"missing.crv";string([base ext])];
isotope = repmat("F18",3,1);
calibration_factor = [NaN;NaN;.4];
pet_start = ["40";"40";"00:00:40"];
sub = ["01";"";""];
manifest = table(filename,isotope,calibration_factor,pet_start,sub);
manifest.Properties.VariableNames{1} = 'FILENAME';
for ext = {'.csv','.tsv','.xlsx'}
    file = fullfile(testCase.TestData.folder,['manifest' ext{1}]);
    if strcmp(ext{1},'.tsv')
        writetable(manifest,file,'FileType','text','Delimiter','\t');
    else
        writetable(manifest,file);
    end
    result = swisstrace_convert_batch(file,'cal_dates',"2025-12-31",'cal_values',.5);
    verifyEqual(testCase, result.status, ["ok";"error";"ok"], strjoin(result.message,newline));
    verifyTrue(testCase, isfile(fullfile(testCase.TestData.folder,'BIDS','sub-01','pet','sub-01_recording-autosampler_blood.tsv')));
    log = readtable(fullfile(testCase.TestData.folder,'calibration_factors.tsv'),'FileType','text','Delimiter','\t');
    verifyEqual(testCase, height(log), 2);
    verifyEqual(testCase, log.gap_days(1), 1);
end
end

function testBatchTableHalfLifeAndEmpty(testCase)
filename = string(exact_file(testCase));
manifest = table(filename,6586.2,.5,'VariableNames',{'filename','half_life','calibration_factor'});
r = swisstrace_convert_batch(manifest,'output_folder',testCase.TestData.folder);
verifyEqual(testCase, r.status, "ok", strjoin(r.message,newline));
r = swisstrace_convert_batch(table(strings(0,1),'VariableNames',{'filename'}),'output_folder',testCase.TestData.folder);
verifyEqual(testCase, height(r), 0);
end

function testRCorrectionParity(testCase)
folder = fullfile(testCase.TestData.root,'tests','generated');
assumeTrue(testCase,isfile(fullfile(folder,'reference.crv')), 'Run generate_r_reference.R for direct R comparisons (CI does this automatically).');
for mode = {'native','auto','framed','nozero','c11'}
    mode = mode{1};
    options = {'pet_start',40.25}; iso = 'F18';
    if strcmp(mode,'auto'), options = {}; end
    if strcmp(mode,'framed'), options = [options,{'frame_scheme',[1 180;10 600]}]; end
    if strcmp(mode,'nozero'), options = [options,{'zero_first_frame',false}]; end
    if strcmp(mode,'c11'), iso = 'C11'; end
    r = swisstrace_correct(fullfile(folder,'reference.crv'),.425,iso,options{:});
    expected = readtable(fullfile(folder,[mode '_tac.tsv']),'FileType','text','Delimiter','\t');
    verifyEqual(testCase, r.tac{:,:}, expected{:,:}, mode, 'AbsTol', 1e-9);
    metadata = readtable(fullfile(folder,[mode '_meta.tsv']),'FileType','text','Delimiter','\t','TextType','string');
    for k = 1:height(metadata)
        value = r.(metadata.name(k));
        if isnan(metadata.value(k))
            verifyTrue(testCase, isnan(value), char(metadata.name(k)));
        else
            verifyEqual(testCase, value, metadata.value(k), char(metadata.name(k)), 'AbsTol', 1e-9);
        end
    end
end
end

function testRExportAndLookupParity(testCase)
folder = fullfile(testCase.TestData.root,'tests','generated');
assumeTrue(testCase,isfile(fullfile(folder,'reference.crv')), 'Run generate_r_reference.R to enable this comparison.');
p = swisstrace_process(fullfile(folder,'reference.crv'),.425,'F18','pet_start',40.25, ...
    'output_folder',testCase.TestData.folder,'sub','01','ses','02');
expected = readmatrix(fullfile(folder,'exports','Corrected_PMOD','reference_corrected.crv'),'FileType','text','NumHeaderLines',1);
actual = readmatrix(p.corrected_crv,'FileType','text','NumHeaderLines',1);
verifyEqual(testCase, actual, expected, 'AbsTol', 1e-9);
expected = readmatrix(fullfile(folder,'exports','BIDS','sub-01','ses-02','pet','sub-01_ses-02_recording-autosampler_blood.tsv'),'FileType','text','NumHeaderLines',1);
actual = readmatrix(p.bids_tsv,'FileType','text','NumHeaderLines',1);
verifyEqual(testCase, actual, expected, 'AbsTol', 1e-9);
a = lookup_calibration(["2025-12-31","2026-01-02"],[.4,.5],'date',["2026-01-01","2026-01-03"]);
b = readtable(fullfile(folder,'lookup.tsv'),'FileType','text','Delimiter','\t');
verifyEqual(testCase, a.calibration_factor, b.calibration_factor);
verifyEqual(testCase, a.gap_days, b.gap_days);
end

function file = exact_file(testCase)
t = (0:300)';
counts = repmat(30,size(t));
post = t >= 40;
counts(post) = counts(post) + 100*exp(-log(2)/6586.2*(t(post)-40));
file = write_raw(testCase.TestData.folder,t,counts);
end

function file = write_raw(folder,t,counts,start)
if nargin < 4, start = datetime(2026,1,1,'TimeZone','UTC'); end
file = fullfile(folder,'raw.crv');
fid = fopen(file,'wt');
cleanup = onCleanup(@() fclose(fid));
for k = 1:numel(t)
    stamp = start + seconds(t(k));
    fprintf(fid,'%d %d %d %d %d %.6f %.17g 550 1000\n', ...
        year(stamp),month(stamp),day(stamp),hour(stamp),minute(stamp),second(stamp),counts(k));
end
end

function touch(file)
fid = fopen(file,'wt'); fclose(fid);
end
