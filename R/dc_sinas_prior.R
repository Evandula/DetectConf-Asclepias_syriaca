#' Build national-level priors from the SInAS database
#'
#' Wraps section 22B of the DetectConf pipeline. Filters the Seebens et
#' al. (2017) SInAS database to records of the focal species marked as
#' introduced and present, then aggregates per country into two
#' rasterised priors:
#' \describe{
#'   \item{\code{national_presence}}{binary 0/1 flag — is the species
#'     known to be introduced and present in this country?}
#'   \item{\code{confirmation_prior}}{continuous 0-1 score combining
#'     introduction status, establishment status, and number of source
#'     datasets, normalised globally.}
#' }
#'
#' These rasters are extracted per \code{cell_id} in
#' \code{dc_assemble_model_frame()} and \code{dc_project_country()},
#' supplying the national-level component of the model's knowledge
#' priors.
#'
#' @param sinas_data A data frame or data.table containing SInAS
#'   records. The full SInAS database (~120 MB) or a pre-subsetted file
#'   (such as the one shipped in \code{inst/extdata/} for
#'   \emph{Asclepias syriaca}) are both accepted. Must contain columns
#'   \code{taxon}, \code{establishmentMeans}, \code{occurrenceStatus},
#'   \code{degreeOfEstablishment}, \code{datasetName}, \code{location},
#'   \code{eventDate}.
#' @param focal_species Scientific name of the focal species with space
#'   separator (e.g. \code{"Asclepias syriaca"}).
#' @param raster_template The classified climate envelope raster (output
#'   of \code{dc_climate_envelope}), used as the rasterisation target.
#' @param world_polygons An \code{sf} or \code{SpatVector} of country
#'   polygons. If \code{NULL}, fetched from
#'   \code{geodata::world(resolution = 1)}. Country names are derived
#'   from \code{GID_0} via \code{countrycode}.
#' @param baseline_prior Numeric baseline value assigned to countries
#'   with no SInAS records. Default \code{0.1}.
#' @param out_dir Optional directory in which to write
#'   \code{sinas_confirmation_prior.tif} and
#'   \code{sinas_national_presence.tif}. If \code{NULL}, rasters are
#'   only returned in memory.
#'
#' @return A named list with two \code{terra::SpatRaster} objects:
#'   \code{confirmation_prior} and \code{national_presence}.
#'
#' @details
#' The \code{confirmation_prior} score is computed per country as
#' \eqn{0.4 \cdot \text{intro\_score} + 0.4 \cdot \text{estab\_score}
#' + 0.2 \cdot \text{evidence\_score}}, where the components are mean
#' introduction-flag, mean established-flag, and \code{log1p} of the
#' number of source datasets. The score is then normalised so the
#' maximum across countries equals 1.
#'
#' Countries not present in SInAS for the focal species receive
#' \code{national_presence = 0} and \code{confirmation_prior =
#' baseline_prior}. This makes the prior a soft floor rather than a
#' hard exclusion — a country with no SInAS record can still be
#' assigned high detection confidence if local effort and corroboration
#' support it.
#'
#' @references
#' Seebens, H., Blackburn, T. M., Dyer, E. E., et al. (2017). No
#' saturation in the accumulation of alien species worldwide.
#' \emph{Nature Communications}, 8, 14435.
#' \doi{10.1038/ncomms14435}
#'
#' @importFrom data.table fread setDT
#' @importFrom terra rasterize vect writeRaster
#' @importFrom sf st_as_sf st_make_valid
#' @importFrom countrycode countrycode
#' @export
dc_sinas_prior <- function(sinas_data,
                           focal_species,
                           raster_template,
                           world_polygons = NULL,
                           baseline_prior = 0.1,
                           out_dir        = NULL) {

  taxon <- establishmentMeans <- occurrenceStatus <- NULL
  degreeOfEstablishment <- datasetName <- location <- eventDate <- NULL
  intro_flag <- estab_flag <- n_sources <- evidence_score <- NULL
  national_presence <- confirmation_prior <- NULL
  intro_score <- estab_score <- year_first_national <- NULL

  if (!data.table::is.data.table(sinas_data))
    sinas_data <- data.table::as.data.table(sinas_data)

  # Treat empty strings and "NULL" as NA
  for (col in names(sinas_data)) {
    set_na <- which(sinas_data[[col]] == "" |
                    sinas_data[[col]] == "NULL")
    if (length(set_na))
      data.table::set(sinas_data, i = set_na, j = col, value = NA)
  }

  # ── Filter to focal species, introduced, present ─────────────
  sinas_spp <- sinas_data[
    taxon == focal_species &
      grepl("introduced", establishmentMeans, ignore.case = TRUE) &
      occurrenceStatus == "present"]

  message("SInAS records for ", focal_species, ": ", nrow(sinas_spp))
  if (nrow(sinas_spp) == 0)
    warning("Species not found in SInAS. All priors will use baseline.")

  if (nrow(sinas_spp) > 0) {
    sinas_spp[, intro_flag := data.table::fifelse(
      grepl("introduced", establishmentMeans, ignore.case = TRUE), 1, 0)]
    sinas_spp[, estab_flag := data.table::fifelse(
      degreeOfEstablishment %in% c("established", "reproducing", "invasive"),
      1, 0)]
    sinas_spp[, n_sources      := lengths(strsplit(datasetName, ";"))]
    sinas_spp[, evidence_score := log1p(n_sources)]

    prior_dt <- sinas_spp[, .(
      national_presence   = 1L,
      intro_score         = mean(intro_flag,     na.rm = TRUE),
      estab_score         = mean(estab_flag,     na.rm = TRUE),
      evidence_score      = mean(evidence_score, na.rm = TRUE),
      year_first_national = {
        y <- suppressWarnings(
          min(as.integer(sub("-.*", "", eventDate)), na.rm = TRUE))
        if (!is.finite(y)) NA_integer_ else as.integer(y)
      }
    ), by = location]

    prior_dt[, confirmation_prior :=
               0.4 * intro_score + 0.4 * estab_score +
               0.2 * evidence_score]
    max_prior <- max(prior_dt$confirmation_prior, na.rm = TRUE)
    if (is.finite(max_prior) && max_prior > 0)
      prior_dt[, confirmation_prior := confirmation_prior / max_prior]
  } else {
    prior_dt <- data.table::data.table(
      location           = character(),
      national_presence  = integer(),
      confirmation_prior = numeric())
  }

  message("Countries with SInAS prior data: ", nrow(prior_dt))

  # ── World polygons ───────────────────────────────────────────
  if (is.null(world_polygons)) {
    world_v <- geodata::world(resolution = 1, path = tempdir())
  } else if (inherits(world_polygons, "sf")) {
    world_v <- terra::vect(world_polygons)
  } else {
    world_v <- world_polygons
  }
  world_sf <- sf::st_as_sf(world_v)
  world_sf$location <- suppressWarnings(
    countrycode::countrycode(world_sf$GID_0, "iso3c", "country.name"))
  world_sf <- sf::st_make_valid(world_sf)

  # ── Merge priors onto world ──────────────────────────────────
  world_sinas <- merge(
    world_sf[, c("location", "geometry")],
    prior_dt[, .(location, national_presence, confirmation_prior)],
    by = "location", all.x = TRUE)

  world_sinas$national_presence[is.na(world_sinas$national_presence)]   <- 0L
  world_sinas$confirmation_prior[is.na(world_sinas$confirmation_prior)] <-
    baseline_prior

  # ── Rasterise ─────────────────────────────────────────────────
  conf_prior_rast   <- terra::rasterize(terra::vect(world_sinas),
                                        raster_template,
                                        field = "confirmation_prior")
  nat_presence_rast <- terra::rasterize(terra::vect(world_sinas),
                                        raster_template,
                                        field = "national_presence")

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(conf_prior_rast,
                       file.path(out_dir, "sinas_confirmation_prior.tif"),
                       overwrite = TRUE)
    terra::writeRaster(nat_presence_rast,
                       file.path(out_dir, "sinas_national_presence.tif"),
                       overwrite = TRUE)
    message("SInAS prior rasters written to: ", out_dir)
  }

  list(confirmation_prior = conf_prior_rast,
       national_presence  = nat_presence_rast)
}


#' Extract raster values by cell ID
#'
#' Convenience wrapper around \code{terra::values} that returns the
#' value of \code{rast_obj} at each \code{cell_id}, preserving order.
#' Used internally by DetectConf to attach raster-based variables to
#' presence, pseudo-absence, and projection data frames.
#'
#' @param rast_obj A \code{terra::SpatRaster} (single layer).
#' @param cell_ids Integer vector of cell indices.
#'
#' @return Numeric vector of length \code{length(cell_ids)}.
#'
#' @importFrom terra values
#' @export
dc_extract_by_cell <- function(rast_obj, cell_ids) {
  as.vector(terra::values(rast_obj))[cell_ids]
}
