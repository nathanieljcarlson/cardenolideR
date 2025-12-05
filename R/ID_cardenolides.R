#' Identify cardenolides from HPLC chromatograms
#'
#' Reads chromatograms, preprocesses them, optionally performs PTW warping,
#' identifies candidate cardenolides based on absorption spectra, and allows
#' manual verification with notes. Exports CSV files with results.
#'
#' @title Identify cardenolides from HPLC chromatograms
#' @description Tools for detecting cardenolides or similar compounds from HPLC data.
#' @param input_dir Directory with raw chromatogram files
#' @param output_dir Directory to save CSV results
#' @param format_in Input chromatogram format (default "chemstation_uv")
#' @param peak_lambda Wavelength to focus on (default 218)
#' @param peak_range Range of absorbance maxima to filter candidate peaks (default 216-222)
#' @param do_ptw Whether to perform PTW retention time warping (default TRUE)
#' @param do_baseline Whether to perform baseline correction (default FALSE)
#' @param do_manual_verification Whether to manually verify peaks (default FALSE)
#' @param minor_peak_max Maximum size of minor peak (outside target lambda) for inclusion (default 0.25)
#' @param skip_auto_screen Skip automatic screening and go straight to manual (default FALSE)
#' @param rt_range_min Minimum retention time to consider (in minutes, default 0)
#' @param rt_range_max Maximum retention time to consider (in minutes, default 40)
#' @param new_ts Retention time sequence for preprocessing
#' @param new_lambdas Wavelength sequence for preprocessing
#' @param n_cores Number of CPU cores for parallel processing (default 4)
#' @return A list containing:
#'   \itemize{
#'     \item all_peaks - data.frame of all detected peaks with sample IDs as row names
#'     \item auto_candidates - data.frame of automatically selected candidate peaks with sample IDs as row names
#'     \item final_peaks - data.frame of peaks identified as cardenolides with sample IDs as row names
#'     \item verified_df - manual verification data (if manual verification used)
#'   }
#' @importFrom plotly %>%
#' @importFrom plotly layout
#' @importFrom utils write.csv
#' @export
ID_cardenolides <- function(
    input_dir = getwd(),
    output_dir = getwd(),
    format_in = "chemstation_uv",
    peak_lambda = 218,
    peak_range = c(216, 222),
    do_ptw = TRUE,
    do_baseline = FALSE,
    do_manual_verification = FALSE,
    minor_peak_max = 0.25,
    skip_auto_screen = FALSE,
    rt_range_min = 0,
    rt_range_max = 40,
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
  x <- chromatographR::read_chroms(input_dir, format_in = format_in, cl = n_cores)

  if(do_baseline){
    message("Applying baseline correction...")
    x <- chromatographR::baseline_correct(x)
  }


  #-------------------------
  # STEP 2: PREPROCESS
  #-------------------------
  message("Preprocessing chromatograms...")
  dat.pr <- chromatographR::preprocess(x, new_ts, new_lambdas)

  #-------------------------
  # STEP 4: PTW WARPING (OPTIONAL)
  #-------------------------
  if(do_ptw){
    message("Building PTW warping models...")
    warping.models <- chromatographR::correct_rt(dat.pr, alg="ptw", what="models", lambdas = peak_lambda, scale = TRUE)
    message("Applying warping models...")
    warp <- chromatographR::correct_rt(dat.pr, models = warping.models)
  } else {
    message("Skipping PTW warp...")
    warp <- dat.pr
    warping.models <- NULL
  }

  #-------------------------
  # STEP 5: READ AND FILTER PEAK LISTS
  #-------------------------
  message("Reading peak reports...")
  pks <- chromConverter::read_peaklist(input_dir)

  #-------------------------
  # STEP 6: CORRECT PEAK RETENTION TIMES (OPTIONAL)
  #-------------------------
  if(do_ptw){
    message("Correcting peak retention times...")
    pks.cor <- chromatographR::correct_peaks(pks, warping.models, chrom_list = dat.pr)
    peaks_for_table <- pks.cor
  } else {
    peaks_for_table <- pks
  }

  #-------------------------
  # STEP 7: BUILD PEAK TABLE
  #-------------------------
  message("Building peak table...")
  pktab <- chromatographR::get_peaktable(peaks_for_table, use.cor = do_ptw)
  pktab$args$chrom_list <- "warp"
  pktab <- chromatographR::filter_peaktable(pktab, lambda = as.character(peak_lambda))

  # Double-check that all peaks are within RT range (safety check)
  rt_values <- as.numeric(pktab[[2]][3, ])
  peaks_in_range <- rt_values >= rt_range_min & rt_values <= rt_range_max

  if(!all(peaks_in_range)) {
    warning("Some peaks are outside the specified RT range. Filtering them out.")
    pktab[[1]] <- pktab[[1]][, peaks_in_range, drop = FALSE]
    pktab[[2]] <- pktab[[2]][, peaks_in_range, drop = FALSE]
    rt_values <- rt_values[peaks_in_range]
  }

  message("Found ", ncol(pktab[[1]]), " peaks within the retention time range")

  # Create data frame with row names for R object
  all_peaks_df <- pktab[[1]]
  colnames(all_peaks_df) <- rt_values

  # Export CSV with SampleID as first column for file export
  write.csv(data.frame(SampleID = rownames(all_peaks_df), all_peaks_df,
                       check.names = FALSE),
            file.path(output_dir, "all_peaks.csv"), row.names = FALSE)

  # Keep all_peaks in memory WITH ROW NAMES (not SampleID column)
  all_peaks <- all_peaks_df  # This has row names, not SampleID column

  #-------------------------
  # STEP 8: AUTO-SCREEN (OPTIONAL)
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

      spectra_raw <- chromatographR::plot_all_spectra(
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
          area = s[maxima_pos]
        )
        in_range <- maxima[maxima$wavelength >= peak_range[1] & maxima$wavelength <= peak_range[2], ]
        if(nrow(in_range) == 0) next

        main_peak_area <- max(in_range$area)
        maxima_filtered <- maxima[maxima$area >= main_peak_area * minor_peak_max, ]
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

  # Create candidate output with row names for R object
  candidate_output_df <- pktab[[1]][, candidate_peaks, drop = FALSE]
  colnames(candidate_output_df) <- rt_values[match(candidate_peaks, colnames(pktab[[1]]))]

  # Export CSV with SampleID as first column for file export
  write.csv(data.frame(SampleID = rownames(candidate_output_df),
                       candidate_output_df, check.names = FALSE),
            file.path(output_dir, "auto_candidates.csv"), row.names = FALSE)

  # Keep in memory WITH ROW NAMES (not SampleID column)
  candidate_output <- candidate_output_df  # This has row names, not SampleID column

  #-------------------------
  # STEP 9: MANUAL VERIFICATION
  #-------------------------
  if(do_manual_verification){
    verified_df <- data.frame(RT = numeric(), Cardenolide = character(),
                              Notes = character(), stringsAsFactors = FALSE)

    for(peak in candidate_peaks){
      rt <- unlist(pktab[[2]][3, names(pktab[[2]]) %in% peak])

      p <- chromatographR::plot_all_spectra(
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
    final_peaks_names <- candidate_peaks[verified_df$Cardenolide == "Y"]
  } else {
    final_peaks_names <- candidate_peaks
  }

  #-------------------------
  # STEP 10: EXPORT FINAL CARDENOLIDE TABLE
  #-------------------------
  n1 <- names(pktab[[1]])[(names(pktab[[1]]) %in% final_peaks_names)]
  n2 <- names(pktab[[2]])[(names(pktab[[2]]) %in% final_peaks_names)]

  # Create final table with row names for R object
  final_tab_df <- pktab[[1]][, n1, drop = FALSE]
  colnames(final_tab_df) <- as.numeric(pktab[[2]][3, n2, drop = FALSE])

  # Export CSV with SampleID as first column for file export
  write.csv(data.frame(SampleID = rownames(final_tab_df),
                       final_tab_df, check.names = FALSE),
            file.path(output_dir, "raw_cards.csv"), row.names = FALSE)

  # Keep in memory WITH ROW NAMES (not SampleID column)
  final_tab <- final_tab_df  # This has row names, not SampleID column

  message("Done! Results saved to: ", output_dir)

  #-------------------------
  # RETURN RESULTS
  #-------------------------
  invisible(list(
    all_peaks = all_peaks,           # Has row names, not SampleID column
    auto_candidates = candidate_output,  # Has row names, not SampleID column
    final_peaks = final_tab,         # Has row names, not SampleID column
    verified_df = if(do_manual_verification) verified_df else NULL
  ))
}
