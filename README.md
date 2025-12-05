# cardenolideR


**An R package for identifying cardenolides from HPLC chromatogram files.**  
=======
cardenolideR leverages [chromatographR](https://ethanbass.github.io/chromatographR/) to read, preprocess, and analyze chromatograms, allowing both automated and manual identification of candidate cardenolides.  

---

## Installation

You can install the development version directly from GitHub:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install cardenolideR from GitHub
devtools::install_github("nathanieljcarlson/cardenolideR", dependencies = TRUE)
Dependencies (installed automatically if missing):
```
## Dependencies (installed automatically if missing)
`chromatographR`

`chromConverter`

`plotly`

`parallel`

`htmlwidgets`

## Example Usage

You can try `ID_cardenolides` with the example data included in this package:

```r
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
## Function: `ID_cardenolides()`

Reads chromatograms, preprocesses, optionally performs PTW retention time warping, screens peaks by absorption spectra, and allows manual verification with notes.

## Arguments
| Argument                 | Type           | Default                   | Description                                                    |
| ------------------------ | -------------- | ------------------------- | -------------------------------------------------------------- |
| `input_dir`              | string         | `getwd()`                 | Directory with raw chromatogram files.                         |
| `output_dir`             | string         | `getwd()`                 | Directory to save CSV results.                                 |
| `format_in`              | string         | `"chemstation_uv"`        | Input chromatogram format.                                     |
| `peak_lambda`            | numeric        | `218`                     | Wavelength to focus on.                                        |
| `peak_range`             | numeric vector | `c(216, 240)`             | Range of absorbance maxima to filter candidate peaks.          |
| `do_ptw`                 | logical        | `TRUE`                    | Whether to perform PTW retention time warping.                 |
| `do_baseline`            | logical        | `FALSE`                   | Whether to apply baseline correction.                          |
| `do_manual_verification` | logical        | `FALSE`                   | Whether to manually verify candidate peaks.                    |
| `skip_auto_screen`       | logical        | `FALSE`                   | Skip automatic screening; all peaks go to manual verification. |
| `new_ts`                 | numeric vector | `seq(.01, 39.95, by=.01)` | Retention time sequence for preprocessing.                     |
| `new_lambdas`            | numeric vector | `seq(200, 318, by=2)`     | Wavelength sequence for preprocessing.                         |
| `n_cores`                | integer        | `4`                       | Number of CPU cores for parallel processing.                   |


## Returns
A list containing:

`all_peaks` – table of all detected peaks with SampleID as the first column.

`auto_candidates` – table of automatically selected candidate peaks.

`final_peaks` – table of peaks identified as cardenolides.

`verified_df` – manual verification data (`if do_manual_verification = TRUE`). Includes RT, Cardenolide (Y/N), and Notes columns.

## Output Files
The function writes several CSV files to output_dir:

| File                        | Description                                                                        |
| --------------------------- | ---------------------------------------------------------------------------------- |
| `all_peaks.csv`             | All detected peaks, rows = samples, columns = peak areas.                          |
| `auto_candidates.csv`       | Candidate peaks selected by automatic screening (or all peaks if skipped).         |
| `verified_cardenolides.csv` | Manual verification results with RT, Y/N, and notes (if manual verification used). |
| `raw_cards.csv`             | Final cardenolide table, only peaks confirmed as cardenolides.                     |


## Citation
If you use cardenolideR in your work, please also cite chromatographR:

Bass, E. (2023). chromatographR: Chromatographic Data Analysis Toolset (version 0.7.3). http://doi.org/10.5281/zenodo.6944334

