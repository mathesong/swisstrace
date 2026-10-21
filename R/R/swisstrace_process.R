# Correct a raw twilite .crv and write the outputs to disk:
#   <output_folder>/Corrected_PMOD/<base>_corrected.crv   (corrected curve)
#   <output_folder>/Plots/<base>_corrected.png            (QC plot)
#   <output_folder>/BIDS/sub-XX[/ses-XX]/pet/...          (BIDS PET blood files)
#
# `output_folder` defaults to the folder containing the source .crv. BIDS output is
# produced only when a subject label is given (`sub`); a `ses` without a `sub` is an
# error. If `bids_dir` (a BIDS dataset root) is given, the blood files go there
# instead: beside the subject's _pet.nii[.gz] image when one is found (taking that
# file's entities), otherwise into bids_dir/sub-XX[/ses-XX]/pet/.
# Requires R/swisstrace_correct.R and R/swisstrace_qc.R to be sourced.

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Correct a twilite .crv and write its outputs (corrected crv, QC plot, BIDS)
#'
#' @param filename Path to the raw twilite `.crv` to correct.
#' @param calibration_factor Per-study calibration multiplier (counts/sec -> kBq/cc).
#' @param isotope Isotope name for the decay correction (REQUIRED — not stored in the
#'   `.crv`). E.g. "F18". Pass `half_life` via `...` instead for an isotope not in the
#'   built-in table.
#' @param ... Further correction arguments passed to [swisstrace_correct()]
#'   (e.g. `pet_start`, `lead`, `half_life`, `frame_scheme`).
#' @param output_folder Root for outputs. Default: the folder containing `filename`.
#'   Outputs are written to subfolders within it.
#' @param sub,ses Optional BIDS subject/session labels (e.g. "01", "01"). If `sub` is
#'   given, BIDS PET blood files are written (under `output_folder/BIDS/...` by default,
#'   or into `bids_dir` if supplied). Supplying `ses` without `sub` is an error. With
#'   neither, no BIDS output is made.
#' @param bids_dir Optional path to an existing BIDS dataset root. When given, the BIDS
#'   blood files are placed there instead of under `output_folder`: if a
#'   `sub-<sub>[_ses-<ses>]..._pet.nii[.gz]` image is found within it, the blood files
#'   are written into that image's folder and named from its entities (so the recording
#'   is associated with that scan); otherwise they go into
#'   `bids_dir/sub-<sub>[/ses-<ses>]/pet/`. Requires `sub`.
#' @param recording BIDS `recording-` label for the blood file (default "autosampler").
#' @param overwrite Overwrite existing files (default TRUE).
#'
#' @return (invisibly) a named list of the paths written.
#' @export
swisstrace_process <- function(filename, calibration_factor, isotope, ...,
                                  output_folder = NULL,
                                  sub = NULL, ses = NULL,
                                  bids_dir = NULL,
                                  recording = "autosampler",
                                  overwrite = TRUE) {
  if (missing(isotope)) {
    stop("`isotope` is required (it is not stored in the .crv): e.g. isotope = \"F18\". ",
         "For an isotope outside the built-in table, pass `half_life` via `...`.")
  }
  assert_raw_crv(filename)   # reject a corrected/processed .crv up front

  # fail fast on an invalid BIDS request, before doing any work
  if (is.null(sub) && !is.null(ses)) {
    stop("`ses` was provided without `sub`. A subject label is required for BIDS ",
         "output; provide `sub`, or omit both to skip BIDS.")
  }
  if (is.null(sub) && !is.null(bids_dir)) {
    stop("`bids_dir` was provided without `sub`. A subject label is needed to locate ",
         "the PET image and name the blood files; provide `sub`.")
  }

  res <- swisstrace_correct(filename, calibration_factor = calibration_factor,
                             isotope = isotope, ...)

  output_folder <- output_folder %||% dirname(filename)
  base <- sub("\\.crv$", "", basename(filename), ignore.case = TRUE)

  ## formatters: time as N.N, value as full-precision float
  fmt_time <- function(t) format(t, nsmall = 1, trim = TRUE, scientific = FALSE)
  fmt_val  <- function(v) {
    s <- vapply(v, function(x) format(x, digits = 17, trim = TRUE, scientific = FALSE),
                character(1))
    s[!grepl("\\.", s)] <- paste0(s[!grepl("\\.", s)], ".0")  # 0 -> "0.0"
    s
  }
  tcol <- fmt_time(res$tac$time)
  vcol <- fmt_val(res$tac$activity)
  written <- list()

  ## --- 1. corrected curve ----------------------------------------------------
  cdir <- file.path(output_folder, "Corrected_PMOD")
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  crv_path <- file.path(cdir, paste0(base, "_corrected.crv"))
  header <- "Corrected_&_calibrated_[kBq/cc]_>___time[seconds]\tvalue[kBq/cc]"
  if (overwrite || !file.exists(crv_path)) {
    writeLines(c(header, paste(tcol, vcol, sep = "\t")), crv_path)
  }
  written$corrected_crv <- crv_path

  ## --- 2. QC plot -----------------------------------------------------------
  pdir <- file.path(output_folder, "Plots")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(pdir, paste0(base, "_corrected.png"))
  if (overwrite || !file.exists(png_path)) {
    grDevices::png(png_path, width = 1000, height = 750, res = 110)
    swisstrace_qc(res)
    grDevices::dev.off()
  }
  written$plot_png <- png_path

  ## --- 3. BIDS PET blood files (only when sub is provided) ------------------
  if (!is.null(sub)) {
    clean <- function(x, entity) gsub("[^A-Za-z0-9]", "",
                                      sub(paste0("^", entity, "-"), "", as.character(x)))
    sub_v <- clean(sub, "sub")
    ses_v <- if (!is.null(ses)) clean(ses, "ses") else NULL
    sub_ses_tree <- function(root)            # root/sub-XX[/ses-XX]/pet
      do.call(file.path, as.list(c(root, paste0("sub-", sub_v),
                                   if (!is.null(ses_v)) paste0("ses-", ses_v), "pet")))
    generic_stem <- paste0("sub-", sub_v,     # sub-XX[_ses-XX]_recording-..._blood
                           if (!is.null(ses_v)) paste0("_ses-", ses_v) else "",
                           "_recording-", recording, "_blood")

    if (is.null(bids_dir)) {
      pet_dir <- sub_ses_tree(file.path(output_folder, "BIDS"))
      stem    <- generic_stem
    } else {
      # Look for the subject's (and session's) PET image inside the dataset and put the
      # blood files beside it, inheriting its entities so the recording is associated
      # with that scan. Match against the file name only, anchored, so e.g. sub-1 does
      # not also catch sub-10.
      name_pat <- paste0("^sub-", sub_v,
                         if (!is.null(ses_v)) paste0("_ses-", ses_v) else "",
                         "(_.*)?_pet\\.nii(\\.gz)?$")
      all_pet <- list.files(bids_dir, pattern = "_pet\\.nii(\\.gz)?$",
                            recursive = TRUE, full.names = TRUE)
      hits <- all_pet[grepl(name_pat, basename(all_pet))]
      if (length(hits)) {
        if (length(hits) > 1L) {
          warning(length(hits), " '_pet.nii[.gz]' files matched sub-", sub_v,
                  if (!is.null(ses_v)) paste0(", ses-", ses_v) else "",
                  " in `bids_dir`; using the first:\n  ", hits[1])
        }
        pet_dir <- dirname(hits[1])
        stem    <- paste0(sub("_pet\\.nii(\\.gz)?$", "", basename(hits[1])),
                          "_recording-", recording, "_blood")
      } else {
        pet_dir <- sub_ses_tree(bids_dir)
        stem    <- generic_stem
      }
    }

    dir.create(pet_dir, recursive = TRUE, showWarnings = FALSE)
    tsv_path  <- file.path(pet_dir, paste0(stem, ".tsv"))
    json_path <- file.path(pet_dir, paste0(stem, ".json"))

    if (overwrite || !file.exists(tsv_path)) {
      writeLines(c("time\twhole_blood_radioactivity", paste(tcol, vcol, sep = "\t")),
                 tsv_path)
    }
    # BIDS PET blood sidecar. Only fields we actually have are included; this is
    # autosampler whole-blood data, so plasma and metabolites are unavailable.
    meta <- list(
      PlasmaAvail         = FALSE,
      WholeBloodAvail     = TRUE,
      MetaboliteAvail     = FALSE,
      DispersionCorrected = FALSE,
      time = list(
        Description = "Time relative to time zero (PET scan start), at the frame midpoint.",
        Units = "s"),
      whole_blood_radioactivity = list(
        Description = paste("Whole blood radioactivity concentration: background-subtracted,",
                            "decay-corrected to time zero, and calibrated."),
        Units = "kBq/mL")
    )
    if (overwrite || !file.exists(json_path)) {
      jsonlite::write_json(meta, json_path, auto_unbox = TRUE, pretty = TRUE)
    }
    written$bids_tsv  <- tsv_path
    written$bids_json <- json_path
  }

  message("Wrote:\n  ", paste(unlist(written), collapse = "\n  "))
  invisible(written)
}
