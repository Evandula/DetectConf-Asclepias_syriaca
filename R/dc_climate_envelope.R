#' Build classified climate envelope raster
#'
#' Wraps sections 1-13 of the DetectConf pipeline. Constructs a classified
#' raster (values 0/1/2/3) describing the relationship between every grid
#' cell on Earth and the species' native climate envelope:
#'
#' \describe{
#'   \item{0}{outside both the min-max envelope and Mahalanobis envelope
#'            \emph{and} outside the native range polygon}
#'   \item{1}{inside the Mahalanobis envelope only}
#'   \item{2}{inside both envelopes}
#'   \item{3}{within the native range polygon}
#' }
#'
#' Cells with value 3 are excluded from pseudo-absence sampling so that
#' detection confidence is only estimated for the introduced range.
#'
#' @param occurrences A data frame of species occurrences with columns
#'   \code{decimalLongitude} and \code{decimalLatitude}. Used to derive
#'   the climate envelope inside the native range.
#' @param native_range_polygon An \code{sf} or \code{SpatVector} polygon
#'   describing the species' native range. The Mahalanobis distance and
#'   min-max ranges are computed from cells inside this polygon.
#' @param worldclim_path Directory in which WorldClim data will be cached
#'   (passed to \code{geodata::worldclim_global}).
#' @param bio_vars Integer vector of WorldClim bioclim variable indices to
#'   consider before VIF filtering. Default is \code{c(5, 6, 10, 11, 13, 14,
#'   16, 17)} (temperature extremes and precipitation extremes).
#' @param vif_threshold Variance inflation factor cutoff for
#'   \code{usdm::vifstep}. Default \code{10}.
#' @param mahal_quantile Quantile of the chi-square distribution used to
#'   threshold Mahalanobis distance. Default \code{0.95}.
#' @param out_file Optional path to write the classified raster as a
#'   GeoTIFF. If \code{NULL}, the raster is only returned in memory.
#'
#' @return A \code{terra::SpatRaster} with integer values in
#'   \{0, 1, 2, 3, NA\}.
#'
#' @details
#' The classified raster is the spatial backbone of the entire DetectConf
#' pipeline. All later functions (\code{dc_extract_effort},
#' \code{dc_generate_pseudo_absences}, \code{dc_project_country}) reference
#' grid cells via \code{cell_id}, which is the cell index in this raster.
#' Save the raster and reuse the same object across the workflow to ensure
#' \code{cell_id} values remain consistent.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' range_pol <- st_read(system.file("extdata",
#'                                  "asclepias_syriaca_range.gpkg",
#'                                  package = "DetectConf"))
#' occ <- read.csv("asclepias_syriaca_gbif.csv", sep = "\t")
#' env <- dc_climate_envelope(occ, range_pol,
#'                            worldclim_path = "data/",
#'                            out_file = "climate_envelope.tif")
#' }
#'
#' @importFrom terra rast values ifel mosaic mask extend rasterize crs
#' @importFrom terra extract vect writeRaster
#' @importFrom geodata worldclim_global
#' @importFrom usdm vifstep
#' @importFrom stats cov mahalanobis qchisq complete.cases
#' @export
dc_climate_envelope <- function(occurrences,
                                native_range_polygon,
                                worldclim_path = "data/",
                                bio_vars       = c(5, 6, 10, 11,
                                                   13, 14, 16, 17),
                                vif_threshold  = 10,
                                mahal_quantile = 0.95,
                                out_file       = NULL) {

  # ── Coerce native range to SpatVector ────────────────────────
  if (inherits(native_range_polygon, "sf")) {
    native_range_polygon <- terra::vect(native_range_polygon)
  }

  # ── Sections 1-2: WorldClim + VIF ────────────────────────────
  bio_all  <- geodata::worldclim_global(var  = "bio",
                                        res  = 5,
                                        path = worldclim_path)
  names(bio_all) <- paste0("BIO", 1:19)
  bio_sel  <- bio_all[[bio_vars]]
  vif      <- usdm::vifstep(bio_sel, th = vif_threshold)
  retained <- setdiff(vif@variables, vif@excluded)
  bio_fin  <- bio_all[[retained]]
  rm(bio_all)
  message("Bio variables retained: ", paste(names(bio_fin), collapse = ", "))

  # ── Section 3: Crop to native range ──────────────────────────
  bio_crop <- terra::mask(bio_fin, native_range_polygon)

  # ── Sections 4-10: Climate envelope ──────────────────────────
  pts_vc <- terra::vect(occurrences,
                        geom = c("decimalLongitude", "decimalLatitude"),
                        crs  = "EPSG:4326")

  sample_env        <- terra::extract(bio_crop, pts_vc)
  sample_env_values <- sample_env[, -1]
  valid_rows        <- stats::complete.cases(sample_env_values)
  sample_env_model  <- sample_env_values[valid_rows, ]

  env_ranges <- apply(sample_env_model, 2, range)
  centroid   <- colMeans(sample_env_model)
  cov_mat    <- stats::cov(sample_env_model)

  all_env   <- as.data.frame(terra::values(bio_fin))
  valid     <- stats::complete.cases(all_env)
  all_valid <- all_env[valid, ]

  mahal_d       <- stats::mahalanobis(all_valid, centroid, cov_mat)
  threshold     <- stats::qchisq(mahal_quantile, df = ncol(sample_env_model))
  similar_cells <- mahal_d < threshold

  analog_rast               <- terra::rast(bio_fin[[1]])
  vals                      <- rep(NA, terra::ncell(analog_rast))
  vals[valid]               <- similar_cells
  terra::values(analog_rast) <- vals

  within_range <- apply(all_valid, 1, function(row)
    all(row >= env_ranges[1, ] & row <= env_ranges[2, ]))

  envelope_rast               <- terra::rast(bio_fin[[1]])
  vals2                       <- rep(NA, terra::ncell(envelope_rast))
  vals2[valid]                <- within_range
  terra::values(envelope_rast) <- vals2

  # ── Section 11: Classified raster (0/1/2/3) ──────────────────
  mean_r <- analog_rast + envelope_rast
  pol_r  <- terra::rasterize(native_range_polygon, mean_r)
  pol_r  <- terra::extend(pol_r, mean_r)
  mean_r <- terra::mosaic(mean_r, pol_r, fun = "mean")
  mean_r <- terra::ifel(!is.na(mean_r) & !(mean_r %in% c(0, 1, 2)),
                        3, mean_r)

  if (!is.null(out_file)) {
    terra::writeRaster(mean_r, out_file, overwrite = TRUE)
    message("Climate envelope written to: ", out_file)
  }

  mean_r
}
