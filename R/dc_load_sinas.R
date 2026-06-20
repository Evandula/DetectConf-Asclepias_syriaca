#' Load the SInAS subset bundled with DetectConf
#'
#' Returns the Asclepias syriaca records extracted from the SInAS
#' alien species database. This subset is provided for the demo only.
#' To apply DetectConf to other species, download the full SInAS
#' database from Zenodo (see references).
#'
#' @details
#' Users of this subset MUST cite the original SInAS publication and
#' the Zenodo deposit, regardless of whether they work with the
#' bundled subset or download the full database:
#'
#' - Seebens, H. et al. (2017). No saturation in the accumulation of
#'   alien species worldwide. Nature Communications 8, 14435.
#' - SInAS v3.1.1: https://doi.org/10.5281/zenodo.17727120
#'
#' @return A data.table containing the SInAS records for
#'   Asclepias syriaca.
#' @export
dc_load_sinas <- function() {
  path <- system.file("extdata", "sinas_asclepias_syriaca_subset.csv",
                       package = "detectconf")
  if (!file.exists(path)) {
    stop("SInAS subset not found. Reinstall the package.")
  }
  data.table::fread(path)
}
