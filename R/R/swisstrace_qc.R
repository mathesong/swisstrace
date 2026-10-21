# QC plot for a twilite correction: raw coincidences (with detected time 0, onset
# and background region) over the corrected/calibrated curve, sharing a time axis
# referenced to time 0. Lets you eyeball the auto-detection per study.

#' QC plot of a twilite correction
#'
#' @param res A `swisstrace_correct()` result list.
#' @param xlim Optional time range (s, relative to time 0). Defaults to the full
#'   baseline through the end of the corrected curve.
#' @param col_raw,col_corr Colours for the raw and corrected traces.
#' @return `res`, invisibly (called for the plot).
#' @export
swisstrace_qc <- function(res, xlim = NULL,
                            col_raw = "grey40", col_corr = "firebrick") {
  stopifnot(is.list(res), !is.null(res$tac), !is.null(res$raw))
  raw <- res$raw
  tac <- res$tac
  raw_rel   <- raw$time - res$t0_seconds                       # time relative to t0
  onset_rel <- if (isTRUE(res$t0_detected)) res$lead else 0    # onset = t0 + lead
  if (is.null(xlim)) xlim <- c(min(raw_rel), max(tac$time, na.rm = TRUE))

  op <- graphics::par(mfrow = c(2, 1), mar = c(2.5, 4, 2.5, 1), oma = c(2.5, 0, 0, 0))
  on.exit(graphics::par(op), add = TRUE)

  ## panel 1 -- raw coincidence rate
  sel <- raw_rel >= xlim[1] & raw_rel <= xlim[2]
  plot(raw_rel[sel], raw$coincidence[sel], type = "p", pch = 16, cex = 0.4,
       col = col_raw, xlim = xlim, xlab = "", ylab = "coincidences (counts/sec)",
       main = sprintf("%s  |  %s  |  %s, calib = %g",
                      basename(res$file), format(res$date), res$isotope,
                      res$calibration_factor))
  usr <- graphics::par("usr")
  graphics::rect(usr[1], usr[3], onset_rel, usr[4],
                 col = grDevices::adjustcolor("steelblue", 0.12), border = NA)
  graphics::abline(h = res$background, col = "steelblue", lty = 2)
  graphics::abline(v = 0, col = "black", lty = 1)              # time 0
  graphics::abline(v = onset_rel, col = "darkgreen", lty = 3)  # detected onset
  graphics::legend("topright", bty = "n", cex = 0.8,
                   legend = c("time 0", "onset (t0 + lead)",
                              sprintf("background = %.1f (n = %d)",
                                      res$background, res$n_background)),
                   col = c("black", "darkgreen", "steelblue"), lty = c(1, 3, 2))

  ## panel 2 -- corrected & calibrated curve
  plot(tac$time, tac$activity, type = "o", pch = 16, cex = 0.4, col = col_corr,
       xlim = xlim, xlab = "", ylab = "activity (kBq/cc)",
       main = "Corrected & calibrated")
  graphics::abline(h = 0, col = "grey70")
  graphics::abline(v = 0, col = "black", lty = 1)
  graphics::mtext("time relative to time 0 (s)", side = 1, outer = TRUE, line = 1)

  invisible(res)
}
