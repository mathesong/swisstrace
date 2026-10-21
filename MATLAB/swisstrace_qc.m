function fig = swisstrace_qc(res, varargin)
%SWISSTRACE_QC Draw raw/corrected traces and return the figure handle.
%   FIG = SWISSTRACE_QC(RES, 'xlim', [], 'col_raw', [.4 .4 .4],
%       'col_corr', [.70 .13 .13], 'visible', 'on')
%   RES is a SWISSTRACE_CORRECT result. Use 'visible','off' for batch export.
p = inputParser;
p.PartialMatching = false;
addParameter(p, 'xlim', []);
addParameter(p, 'col_raw', [.4 .4 .4]);
addParameter(p, 'col_corr', [.70 .13 .13]);
addParameter(p, 'visible', 'on');
parse(p, varargin{:});
o = p.Results;
raw_rel = res.raw.time - res.t0_seconds;
limits = o.xlim;
if isempty(limits), limits = [min(raw_rel), max(res.tac.time)]; end
if limits(1) == limits(2), limits = limits + [-.5 .5]; end
validateattributes(limits, {'numeric'}, {'real','finite','numel',2,'increasing'});
onset = 0;
if res.t0_detected, onset = res.lead; end
fig = figure('Visible', o.visible, 'Position', [100 100 1000 750], 'Color', 'w');
try
    ax = subplot(2, 1, 1, 'Parent', fig);
    plot(ax, raw_rel, res.raw.coincidence, '.', 'Color', o.col_raw, 'MarkerSize', 4);
    hold(ax, 'on');
    xlim(ax, limits);
    if isfinite(onset)
        yy = ylim(ax);
        patch(ax, [limits(1) onset onset limits(1)], [yy(1) yy(1) yy(2) yy(2)], ...
            'b', 'FaceAlpha', .08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        xline(ax, onset, ':', 'Color', [0 .4 0], 'DisplayName', 'onset (t0 + lead)');
    end
    xline(ax, 0, 'k-', 'DisplayName', 'time 0');
    yline(ax, res.background, '--', 'Color', [.27 .51 .71], ...
        'DisplayName', sprintf('background = %.1f (n = %d)', res.background, res.n_background));
    ylabel(ax, 'coincidences (counts/sec)');
    [~, base, ext] = fileparts(res.file);
    title(ax, sprintf('%s%s | %s | %s, calib = %g', base, ext, ...
        char(string(res.date, 'yyyy-MM-dd')), char(string(res.isotope)), res.calibration_factor), 'Interpreter', 'none');
    legend(ax, 'show', 'Location', 'northeast');
    ax2 = subplot(2, 1, 2, 'Parent', fig);
    plot(ax2, res.tac.time, res.tac.activity, '.-', 'Color', o.col_corr, 'MarkerSize', 4);
    xlim(ax2, limits);
    yline(ax2, 0, '-', 'Color', [.7 .7 .7]);
    xline(ax2, 0, 'k-');
    title(ax2, 'Corrected & calibrated');
    ylabel(ax2, 'activity (kBq/cc)');
    xlabel(ax2, 'time relative to time 0 (s)');
catch err
    close(fig);
    rethrow(err);
end
end
