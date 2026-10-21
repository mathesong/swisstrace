# Generate deterministic raw data and expected outputs with the actual R code.
# Run from the repository root: Rscript MATLAB/tests/generate_r_reference.R
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath('.')
for (f in list.files(file.path(root, 'R', 'R'), pattern = '[.]R$', full.names = TRUE)) source(f)
folder <- file.path(root, 'MATLAB', 'tests', 'generated')
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
t <- seq(0, 300, by = 0.5)
coinc <- 30 + 2 * sin(t)
post <- t >= 60
coinc[post] <- coinc[post] + 150 * (1 - exp(-(t[post] - 60)/2)) * exp(-(t[post] - 60)/300)
stamps <- as.POSIXct('2026-01-01 12:00:00', tz = 'UTC') + t
parts <- as.POSIXlt(stamps, tz = 'UTC')
rows <- sprintf('%d %d %d %d %d %.6f %.17g 550 1000', parts$year+1900, parts$mon+1,
                parts$mday, parts$hour, parts$min, parts$sec, coinc)
file <- file.path(folder, 'reference.crv')
writeLines(rows, file)
write_tsv <- function(x, path) write.table(x, path, sep='\t', quote=FALSE, row.names=FALSE, na='NaN')
for (mode in c('native', 'auto', 'framed', 'nozero', 'c11')) {
    a <- list(file=file, calibration_factor=.425, isotope='F18')
    if (mode != 'auto') a$pet_start <- 40.25
    if (mode == 'framed') a$frame_scheme <- data.frame(width=c(1,10), end=c(180,600))
    if (mode == 'nozero') a$zero_first_frame <- FALSE
    if (mode == 'c11') a$isotope <- 'C11'
    res <- do.call(swisstrace_correct, a)
    write_tsv(res$tac, file.path(folder, paste0(mode, '_tac.tsv')))
    metadata <- unlist(res[c('background','t0_seconds','n_background','half_life','lambda','lead','n_raw')])
    write_tsv(data.frame(name=names(metadata), value=unname(metadata)), file.path(folder, paste0(mode, '_meta.tsv')))
}
lookup <- lookup_calibration(c('2025-12-31','2026-01-02'), c(.4,.5), date=c('2026-01-01','2026-01-03'))
write_tsv(lookup, file.path(folder, 'lookup.tsv'))
swisstrace_process(file, .425, 'F18', pet_start=40.25,
                   output_folder=file.path(folder, 'exports'), sub='01', ses='02')
