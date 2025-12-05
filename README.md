# cardenolideR

An R package for identifying cardenolides from HPLC chromatogram files.
cardenolideR leverages chromatographR to read, preprocess, and analyze
chromatograms, allowing both automated and manual identification of
candidate cardenolides.

## Installation

You can install the development version directly from GitHub:

``` r
# Install devtools if you don't have it
install.packages("devtools")

# Install cardenolideR from GitHub
devtools::install_github("nathanieljcarlson/cardenolideR", dependencies = TRUE)
```

## Dependencies (installed automatically if missing)

`chromatographR`

`chromConverter`

`plotly`

`parallel`

`htmlwidgets`

##Example Usage You can try ID_cardenolides with the example data
included in this package:

``` r
# Load the package
library(cardenolideR)

# Find the path to the example data
example_dir <- system.file("extdata", "example_chroms", package = "cardenolideR")

# Create a temporary output folder
output_dir <- tempdir()

# Run the cardenolide identification
results <- ID_cardenolides(
  input_dir = example_dir,
  output_dir = output_dir,
  do_ptw = TRUE,
  do_manual_verification = TRUE
)

# Check the results
head(results$all_peaks)
head(results$auto_candidates)
head(results$final_peaks)
```

## Function: ID_cardenolides()

Reads chromatograms, preprocesses them, optionally performs PTW
retention time warping, identifies candidate cardenolides based on
absorption spectra, and allows manual verification with notes. Exports
CSV files and interactive HTML absorbance spectra.

## Arguments

| Argument | Type | Default | Description |
|:---|:---|:---|:---|
| `input_dir` | string | `getwd()` | Directory with raw chromatogram files |
| `output_dir` | string | `getwd()` | Directory to save CSV results and HTML spectra |
| `format_in` | string | `"chemstation_uv"` | Input chromatogram format |
| `peak_lambda` | numeric | `218` | Wavelength to focus on for analysis |
| `peak_range` | numeric vector | `c(216, 222)` | Range of absorbance maxima to filter candidate peaks (nm) |
| `do_ptw` | logical | `TRUE` | Whether to perform PTW retention time warping |
| `do_baseline` | logical | `FALSE` | Whether to apply baseline correction |
| `do_manual_verification` | logical | `FALSE` | Whether to manually verify candidate peaks |
| `minor_peak_max` | numeric | `0.25` | Maximum size of minor peak (outside target lambda) for inclusion |
| `skip_auto_screen` | logical | `FALSE` | Skip automatic screening; all peaks go to manual verification |
| `rt_range_min` | numeric | `0` | Minimum retention time to consider (minutes) |
| `rt_range_max` | numeric | `40` | Maximum retention time to consider (minutes) |
| `new_ts` | numeric vector | `seq(.01, 39.95, by=.01)` | Retention time sequence for preprocessing |
| `new_lambdas` | numeric vector | `seq(200, 318, by=2)` | Wavelength sequence for preprocessing |
| `n_cores` | integer | `4` | Number of CPU cores for parallel processing |

## Returns 

A list containing:

`all_peaks` – data.frame of all detected peaks with sample IDs as row
names.

`auto_candidates` – data.frame of automatically selected candidate peaks
with sample IDs as row names.

`final_peaks` – data.frame of peaks identified as cardenolides with
sample IDs as row names.

`verified_df` – manual verification data (only if do_manual_verification
= TRUE). Includes RT, Cardenolide (Y/N), and Notes columns.

`spectra_dir` – directory path containing exported HTML absorbance
spectra.

## Output Files The function writes several files to output_dir:

File Description:

all_peaks.csv - Contains all detected peaks, with rows representing
samples and columns representing peak areas (retention times as column
headers).

auto_candidates.csv - Contains candidate peaks selected by automatic
screening (or all peaks if skipped).

verified_cardenolides.csv - Contains manual verification results with
RT, Y/N, and notes (if manual verification used).

raw_cards.csv - Contains the final cardenolide table, only including
peaks confirmed as cardenolides.

The function also creates a subdirectory called absorbance_spectra/
containing interactive HTML absorbance spectra for all detected peaks.
Each spectrum is named spectrum\_[peak_name]\_RT[time].html.

## Citation 

If you use cardenolideR in your work, please also cite chromatographR:

Bass, E. (2023). chromatographR: Chromatographic Data Analysis Toolset
(version 0.7.3). <http://doi.org/10.5281/zenodo.6944334>
