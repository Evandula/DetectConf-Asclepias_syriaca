#' Download GBIF occurrence data for the demo
#'
#' Downloads the Asclepias syriaca occurrence records used in the
#' DetectConf demo from GBIF. The file is cached locally after the
#' first download and reused on subsequent calls. No GBIF credentials
#' are required.
#'
#' GBIF download DOI: https://doi.org/10.15468/dl.pg5rph
#'
#' @param destdir Directory in which to cache the downloaded file.
#'   Defaults to a persistent user data location.
#' @param overwrite Logical. If TRUE, re-downloads even if a cached
#'   copy exists.
#'
#' @return The file path to the unpacked CSV.
#' @export
dc_download_demo_gbif <- function(destdir = tools::R_user_dir("detectconf", "data"),
                                  overwrite = FALSE) {
  
  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
  
  url      <- "https://api.gbif.org/v1/occurrence/download/request/0071650-260519110011954.zip"
  zip_path <- file.path(destdir, "asclepias_syriaca_gbif.zip")
  csv_path <- file.path(destdir, "0071650-260519110011954.csv")
  
  if (!file.exists(csv_path) || overwrite) {
    message("Downloading GBIF occurrence data (~13 MB) from doi:10.15468/dl.pg5rph...")
    utils::download.file(url, zip_path, mode = "wb")
    utils::unzip(zip_path, exdir = destdir)
    file.remove(zip_path)
  } else {
    message("Using cached GBIF data at: ", csv_path)
  }
  
  return(csv_path)
}