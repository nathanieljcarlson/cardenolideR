#' Identify cardenolides from LC-MS chromatograms
#'
#' Reads chromatograms, preprocesses them, optionally performs PTW warping,
#' identifies candidate cardenolides based on absorption spectra, and allows
#' manual verification with notes. Exports CSV files with results.
#'
#' @title Identify cardenolides from LC-MS chromatograms
#' @description Tools for detecting cardenolides or similar compounds from LC-MS data.
#' @param input_dir Directory with raw chromatogram files
#' @param output_dir Directory to save CSV results
#' @param format_in Input chromatogram format (default "chemstation_uv")
#' @param peak_lambda Wavelength to focus on (default 218)
#' @param peak_range Range of absorbance maxima to filter candidate peaks (default 216-240)
#' @param do_ptw Whether to perform PTW retention time warping (default TRUE)
#' @param do_baseline Whether to perform baseline correction (default FALSE)
#' @param do_manual_verification Whether to manually verify peaks (default FALSE)
#' @param skip_auto_screen Skip automatic screening and go straight to manual (default FALSE)
#' @param new_ts Retention time sequence for preprocessing
#' @param new_lambdas Wavelength sequence for preprocessing
#' @param n_cores Number of CPU cores for parallel processing (default 4)
#' @return A list containing:
#'   \itemize{
#'     \item all_peaks - table of all detected peaks
#'     \item auto_candidates - table of automatically selected candidate peaks
#'     \item final_peaks - table of peaks identified as cardenolides
#'     \item verified_df - manual verification data (if manual verification used)
#'   }
#' @export
ID_cardenolides <- function(
    input_dir = getwd(),
    output_dir = getwd(),
    format_in = "chemstation_uv",
    peak_lambda = 218,
    peak_range = c(216, 240),
    do_ptw = TRUE,
    do_baseline = FALSE,
    do_manual_verification = FALSE,
    skip_auto_screen = FALSE,
    new_ts = seq(.01, 39.95, by=.01),
    new_lambdas = seq(200, 318, by=2),
    n_cores = 4
){

  #-------------------------
  # CREATE OUTPUT DIRECTORY
  #-------------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  #-------------------------
  # STEP 1: READ CHROMATOGRAMS
  #-------------------------
  message("Reading chromatograms from: ", input_dir)
  x <- read_chroms(input_dir, format_in = format_in, cl = n_cores)

  if(do_baseline){
    message("Applying baseline correction...")
    x <- baseline_correct(x)
  }

  #-------------------------
  # STEP 2: PREPROCESS
  #-------------------------
  message("Preprocessing chromatograms...")
  dat.pr <- preprocess(x, new_ts, new_lambdas)

  #-------------------------
  # STEP 3: PTW WARPING (OPTIONAL)
  #-------------------------
  if(do_ptw){
    message("Building PTW warping models...")
    warping.models <- correct_rt(dat.pr, alg="ptw", what="models", lambdas = peak_lambda, scale = TRUE)
    message("Applying warping models...")
    warp <- correct_rt(dat.pr, models = warping.models)
  } else {
    message("Skipping PTW warp...")
    warp <- dat.pr
    warping.models <- NULL
  }

  #-------------------------
  # STEP 4: READ PEAK LISTS
  #-------------------------
  message("Reading peak reports...")
  pks <- read_peaklist(input_dir)

  if(do_ptw){
    message("Correcting peak retention times...")
    pks.cor <- correct_peaks(pks, warping.models, chrom_list = dat.pr)
    peaks_for_table <- pks.cor
  } else {
    peaks_for_table <- pks
  }

  #-------------------------
  # STEP 5: BUILD PEAK TABLE
  #-------------------------
  message("Building peak table...")
  pktab <- get_peaktable(peaks_for_table, use.cor = do_ptw)
  pktab$args$chrom_list <- "warp"
  pktab <- filter_peaktable(pktab, lambda = as.character(peak_lambda))

  # Export all peaks
  all_peaks <- data.frame(SampleID = rownames(pktab[[1]]), pktab[[1]], check.names = FALSE)
  colnames(all_peaks)[-1] <- as.numeric(pktab[[2]][3,])
  write.csv(all_peaks, file.path(output_dir, "all_peaks.csv"), row.names = FALSE)

  #-------------------------
  # STEP 6: AUTO-SCREEN (OPTIONAL)
  #-------------------------
  if(skip_auto_screen){
    message("Skipping automatic candidate selection. All peaks will go to manual verification.")
    candidate_peaks <- colnames(pktab[[1]])
  } else {
    message("Auto-screening spectra for peaks between ", peak_range[1], "-", peak_range[2], " nm...")
    num_peaks <- ncol(pktab[[1]])
    candidate_peaks <- vector("character", num_peaks)

    for(c in seq_len(num_peaks)){
      peak_name <- colnames(pktab[[1]])[c]

      spectra_raw <- plot_all_spectra(
        peak_name,
        peak_table = pktab,
        chrom_list = warp,
        export_spectrum = TRUE,
        scale_spectrum = FALSE,
        verbose = FALSE,
        plot_spectrum = FALSE
      )

      # Ensure spectra is matrix
      if (is.numeric(spectra_raw)) {
        spectra <- matrix(spectra_raw, ncol = 1)
      } else if (is.list(spectra_raw)) {
        spectra <- do.call(cbind, spectra_raw)
      } else {
        stop("Unexpected output from plot_all_spectra()")
      }

      rownames(spectra) <- new_lambdas
      n_samples <- ncol(spectra)
      keep_peak <- FALSE

      for(m in seq_len(n_samples)){
        s <- spectra[, m]
        maxima_pos <- which(diff(sign(diff(s))) == -2) + 1
        maxima <- data.frame(
          wavelength = as.numeric(rownames(spectra)[maxima_pos]),
          height = s[maxima_pos]
        )
        in_range <- maxima[maxima$wavelength >= peak_range[1] & maxima$wavelength <= peak_range[2], ]
        if(nrow(in_range) == 0) next

        main_peak_height <- max(in_range$height)
        maxima_filtered <- maxima[maxima$height >= main_peak_height / 4, ]
        wavelengths <- maxima_filtered$wavelength

        if(length(wavelengths) == 1 &&
           wavelengths >= peak_range[1] &&
           wavelengths <= peak_range[2]){
          keep_peak <- TRUE
        }
      }
      candidate_peaks[c] <- ifelse(keep_peak, peak_name, NA)
    }
    candidate_peaks <- candidate_peaks[!is.na(candidate_peaks)]
  }

  candidate_output <- data.frame(SampleID = rownames(pktab[[1]]), pktab[[1]][, candidate_peaks, drop = FALSE], check.names = FALSE)
  colnames(candidate_output)[-1] <- unlist(pktab[[2]][3, candidate_peaks])
  write.csv(candidate_output, file.path(output_dir, "auto_candidates.csv"), row.names = FALSE)

  #-------------------------
  # STEP 7: MANUAL VERIFICATION
  #-------------------------
  if(do_manual_verification){
    verified_df <- data.frame(RT = numeric(), Cardenolide = character(), Notes = character(), stringsAsFactors = FALSE)

    for(peak in candidate_peaks){
      rt <- unlist(pktab[[2]][3, names(pktab[[2]]) %in% peak])

      p <- plot_all_spectra(
        peak,
        peak_table = pktab,
        chrom_list = warp,
        export_spectrum = FALSE,
        scale_spectrum = FALSE,
        verbose = FALSE,
        plot_spectrum = TRUE,
        engine = "plotly"
      )
      p <- p %>% layout(title = paste0("RT: ", rt, "; Peak: ", peak))
      print(p)

      ans <- readline(prompt = paste0("Is ", peak, " a cardenolide? (Y/N): "))
      ans <- toupper(ans)
      note <- readline(prompt = "Notes? (optional, hit Enter to skip): ")

      verified_df <- rbind(verified_df, data.frame(
        RT = rt,
        Cardenolide = ifelse(ans %in% c("Y","YES"), "Y", "N"),
        Notes = note,
        stringsAsFactors = FALSE
      ))
    }

    write.csv(verified_df, file.path(output_dir, "verified_cardenolides.csv"), row.names = FALSE)
    final_peaks <- candidate_peaks[verified_df$Cardenolide == "Y"]
  } else {
    final_peaks <- candidate_peaks
  }

  #-------------------------
  # STEP 8: EXPORT FINAL CARDENOLIDE TABLE
  #-------------------------
  n1 <- names(pktab[[1]])[(names(pktab[[1]]) %in% final_peaks)]
  n2 <- names(pktab[[2]])[(names(pktab[[2]]) %in% final_peaks)]

  final_tab <- data.frame(SampleID = rownames(pktab[[1]]), pktab[[1]][, n1, drop = FALSE], check.names = FALSE)
  colnames(final_tab)[-1] <- as.numeric(pktab[[2]][3, n2, drop = FALSE])
  write.csv(final_tab, file.path(output_dir, "raw_cards.csv"), row.names = FALSE)

  message("Done! Results saved to: ", output_dir)

  #-------------------------
  # RETURN RESULTS
  #-------------------------
  invisible(list(
    all_peaks = all_peaks,
    auto_candidates = candidate_output,
    final_peaks = final_tab,
    verified_df = if(do_manual_verification) verified_df else NULL
  ))
}

#====================================================================
# USAGE EXAMPLE
#====================================================================
# source("ID_cardenolides.R")
# ID_cardenolides(
#   input_dir = "~/Desktop/2 species total cards/Block 1/K_hel",
#   output_dir = "~/Desktop/processed_cards",
#   do_ptw = TRUE,
#   do_manual_verification = TRUE
# )
